"""
SyGMa comes currently with two rulesets:

phase1
    Phase 1 metabolism rules include mainly different types of oxidation, hydrolysis, reduction
    and condensation reactions

phase2
    Phase 2 metabolism rules include severaly conjugation reaction,
    i.e. with glucuronyl, sulfate, methyl and acetyl
"""
import pkg_resources

ruleset = {
    "phase1": pkg_resources.resource_filename('sygma', "rules/phase1.txt"),
    "phase2": pkg_resources.resource_filename('sygma', "rules/phase2.txt"),
}
