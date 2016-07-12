import pkg_resources

ruleset = {
    "phase1": pkg_resources.resource_filename('sygma', "rules/phase1.txt"),
    "phase2": pkg_resources.resource_filename('sygma', "rules/phase2.txt"),
}
