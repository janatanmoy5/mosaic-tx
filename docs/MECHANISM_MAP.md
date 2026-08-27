# Mechanism map: from donor-specific antibody to graft-injury trajectory

## The chain under test

```text
Donor-specific antibody (DSA)
          ↓
Antibody engagement of donor graft endothelium
          ↓
FCGR3A/CD16-positive NK-cell and monocyte activation
          ↓
Cytotoxicity plus IFN-gamma/chemokine response
          ↓
Microvascular endothelial injury
          ↓
Persistent versus treatment-responsive graft injury
```

This is a sequence of testable propositions, not an assumed causal conclusion. Public expression data
cannot directly observe antibody binding, Fc-receptor cross-linking or cell killing. Those arrows are
supported only when clinical DSA metadata, cell-resolved expression, histology, longitudinal response
and published functional evidence agree.

## Operationalization of every node and arrow

| Node or transition | Observable evidence | Primary resource | Falsification or limitation |
|---|---|---|---|
| DSA exposure | contemporaneous DSA-positive versus DSA-negative label | GSE50084 | DSA titre, subclass and complement-binding information may be unavailable |
| DSA to rejection | DSA+/ABMR versus DSA+/no rejection | GSE50084 | a null result rejects a simple “DSA is sufficient” model |
| Endothelial engagement | DSA context plus glomerular/peritubular endothelial activation | E-MTAB-12051; bulk endothelial scores | RNA does not measure antibody physically bound to endothelium |
| FCGR3A-positive activation | FCGR3A/Fc-signaling program in separately annotated NK and monocyte pseudobulk | E-MTAB-12051; GSE216203 sensitivity | bulk FCGR3A cannot identify lineage; blood NK is not graft NK |
| Cytotoxicity and IFN response | cytotoxic and IFNG/chemokine module effects and GSEA | GSE50084; GSE36059; cell pseudobulk | programs also occur in TCMR and viral infection; specificity must be tested |
| Microvascular injury | endothelial activation/injury module, Banff MVI when available | GSE36059 and cell-resolved endothelial states | expression is not histological injury unless linked to biopsy scores |
| Treatment response | placebo-adjusted patient-blocked change at W24 and W52 | GSE275824 | association with therapy does not identify which cell was directly targeted |
| Persistence | immune module falls while endothelial/injury module remains | GSE275824 | must report confidence intervals; absence of significance is not proof of persistence |

## Prespecified tests

1. **DSA is not sufficient:** compare DSA+/rejection with DSA+/no rejection, retaining a null result.
2. **Fc-effector localization:** score FCGR3A signaling separately in NK and monocyte patient-level
   pseudobulk; adjust interpretation using T-cell and monocyte negative-control modules.
3. **ABMR specificity:** require comparison with both non-rejection and TCMR in GSE36059.
4. **Coupling:** estimate patient/cohort-level association between Fc-effector, cytotoxic/IFNG and
   endothelial-injury scores. Correlation supports coupling but not causal direction.
5. **Temporal uncoupling:** test treatment-by-time interactions for each module. Then test the
   prespecified contrast between immune-module change and endothelial-module change.
6. **Confounding:** repeat analyses excluding or stratifying CMV/BK infection where annotations exist.

## Evidence states

Each arrow is reported as **supported**, **not supported**, or **not testable with available data**.
“Not significant” and “equivalent” are not interchangeable. Formal equivalence bounds must be declared
before claiming absence of a biologically meaningful treatment effect.

## Biomarker and therapeutic consequences

Only genes or modules surviving the complete evidence chain enter biomarker modeling. A biomarker may
combine an active immune component with a persistent endothelial component, but must be validated at
the patient level and tested against TCMR. Drug reversal prioritizes compounds capable of opposing the
supported disease direction. QSAR begins only if that prioritization yields a single tractable target
with adequate comparable activity data.
