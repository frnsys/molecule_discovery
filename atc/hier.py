from .util import concat, split_levels, traverse_to_depth


class ATCHierarchy:
    def __init__(self, atcs):
        # Construct a hierarchy tree
        # for the ATC codes
        self._atcs = atcs
        self._hier = {}
        for cid, codes in atcs.items():
            for code in codes:
                levels = split_levels(code)
                grp = self._hier
                last = levels[-1]
                for l in levels[:-1]:
                    if l not in grp:
                        grp[l] = {}
                    grp = grp[l]
                if last not in grp:
                    grp[last] = set()
                grp[last].add(cid)

    def codes_for_level(self, level):
        return dict(concat(traverse_to_depth(k, children, level)
                           for k, children in self._hier.items()))
