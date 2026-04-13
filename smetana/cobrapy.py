from reframed.core.model import Compartment, Metabolite, ReactionType
from reframed.core.cbmodel import CBModel, CBReaction
from reframed.core.transformation import clean_bounds


def _from_cobrapy_with_prefixes(cb_model):
    """Convert a COBRApy model to a reframed CBModel with R_/M_ ID prefixes.

    Reframed's SBML loader (fbc2 flavor) adds 'R_' to reaction IDs and 'M_'
    to metabolite IDs. COBRApy does not. Since SMETANA's entire codebase
    assumes these prefixes in hardcoded string patterns, we apply them here
    so that converted models are namespace-compatible with SBML-loaded ones.
    """
    try:
        import cobra as cb
    except ImportError:
        raise RuntimeError("COBRApy is not installed.")

    model = CBModel(cb_model.id)

    for c_id, name in cb_model.compartments.items():
        comp = Compartment(c_id, name)
        model.add_compartment(comp)

    for cb_met in cb_model.metabolites:
        m_id = 'M_' + cb_met.id if not cb_met.id.startswith('M_') else cb_met.id
        met = Metabolite(m_id, cb_met.name, cb_met.compartment)
        model.add_metabolite(met)

    for cb_rxn in cb_model.reactions:
        r_id = 'R_' + cb_rxn.id if not cb_rxn.id.startswith('R_') else cb_rxn.id
        stoichiometry = {}
        for cb_met, val in cb_rxn.metabolites.items():
            m_id = 'M_' + cb_met.id if not cb_met.id.startswith('M_') else cb_met.id
            stoichiometry[m_id] = val

        # GPR parsing is skipped: reframed's parse_gpr_rule can hang on complex
        # ModelSEED rules, and SMETANA does not use gene associations.
        has_gpr = bool(cb_rxn.gene_reaction_rule)
        rxn = CBReaction(r_id, name=cb_rxn.name, reversible=(cb_rxn.lower_bound < 0),
                         stoichiometry=stoichiometry, lb=cb_rxn.lower_bound, ub=cb_rxn.upper_bound,
                         objective=cb_rxn.objective_coefficient)
        rxn._has_gpr = has_gpr
        model.add_reaction(rxn)

    return model


def _detect_reaction_types(model):
    """Detect and assign reaction types on a CBModel converted from COBRApy.

    Reframed's from_cobrapy() leaves all reactions as ReactionType.OTHER.
    This function classifies them using the same logic as reframed's SBML
    loader (unbalanced-metabolite detection), without requiring the libsbml
    object.
    """
    for r_id, rxn in model.reactions.items():
        substrates = [m_id for m_id, coeff in rxn.stoichiometry.items() if coeff < 0]
        products = [m_id for m_id, coeff in rxn.stoichiometry.items() if coeff > 0]

        # Exchange: unbalanced reaction in the external compartment
        if len(substrates) == 0 or len(products) == 0:
            compartments = {model.metabolites[m_id].compartment
                            for m_id in rxn.stoichiometry}
            if any(model.compartments[c].external for c in compartments):
                rxn.reaction_type = ReactionType.EXCHANGE
            else:
                rxn.reaction_type = ReactionType.SINK
            continue

        # Transport: spans multiple compartments
        compartments = {model.metabolites[m_id].compartment
                        for m_id in rxn.stoichiometry}
        if len(compartments) > 1:
            rxn.reaction_type = ReactionType.TRANSPORT
            continue

        # Enzymatic: has gene associations
        if getattr(rxn, '_has_gpr', False) or rxn.gpr is not None:
            rxn.reaction_type = ReactionType.ENZYMATIC
            continue

        rxn.reaction_type = ReactionType.OTHER


def _detect_external_compartment(model):
    """Mark the extracellular compartment as external.

    Uses the COBRApy convention: the compartment named 'e' or 'extracellular',
    or failing that, the compartment that participates in the most unbalanced
    (boundary) reactions.
    """
    # Try common COBRApy extracellular compartment IDs
    for c_id in ('e', 'e0', 'extracellular'):
        if c_id in model.compartments:
            model.compartments[c_id].external = True
            return

    # Fallback: find the compartment with the most boundary reactions
    boundary_compartments = []
    for rxn in model.reactions.values():
        substrates = [m for m, c in rxn.stoichiometry.items() if c < 0]
        products = [m for m, c in rxn.stoichiometry.items() if c > 0]
        if len(substrates) == 0 or len(products) == 0:
            for m_id in rxn.stoichiometry:
                boundary_compartments.append(model.metabolites[m_id].compartment)

    if boundary_compartments:
        ext_comp = max(set(boundary_compartments), key=boundary_compartments.count)
        model.compartments[ext_comp].external = True


def convert_cobrapy_model(cobra_model):
    """Convert a single COBRApy model to a reframed CBModel ready for SMETANA.

    Performs the following beyond what reframed's from_cobrapy() does:
      - Adds R_/M_ prefixes to reaction/metabolite IDs to match the BiGG
        namespace that SMETANA's internals expect
      - Detects and marks the extracellular compartment
      - Classifies reaction types (exchange, sink, transport, enzymatic)
      - Zeros the ATP maintenance reaction lower bound

    Args:
        cobra_model: a cobra.Model object

    Returns:
        CBModel: a reframed model compatible with SMETANA
    """
    model = _from_cobrapy_with_prefixes(cobra_model)

    _detect_external_compartment(model)
    _detect_reaction_types(model)
    clean_bounds(model)

    if 'R_ATPM' in model.reactions:
        model.reactions.R_ATPM.lb = 0

    return model


def convert_cobrapy_models(cobra_models):
    """Convert a list of COBRApy models to reframed CBModels.

    Args:
        cobra_models: list of cobra.Model objects

    Returns:
        list of CBModel objects ready for SMETANA
    """
    return [convert_cobrapy_model(m) for m in cobra_models]


class CobraModelCache:
    """A ModelCache-compatible wrapper for pre-converted COBRApy models.

    Implements the same interface as reframed's ModelCache (get_ids, get_model)
    so it can be used as a drop-in replacement in the SMETANA pipeline.
    """

    def __init__(self, cobra_models):
        """
        Args:
            cobra_models: list of cobra.Model objects
        """
        self._models = {}
        for cobra_model in cobra_models:
            model = convert_cobrapy_model(cobra_model)
            self._models[model.id] = model

    def get_ids(self):
        return list(self._models.keys())

    def get_model(self, model_id, reset_id=False):
        if model_id not in self._models:
            raise RuntimeError("Model not in list: " + model_id)
        model = self._models[model_id]
        if reset_id:
            model.id = model_id
        return model
