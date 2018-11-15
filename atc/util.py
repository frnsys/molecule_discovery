def concat(lists):
    """Concat a group of lists"""
    return sum(lists, [])


def get_leaves(node):
    """Get descendant leaves for a node"""
    if isinstance(node, dict):
        return concat(map(get_leaves, node.values()))
    else:
        return list(node)


def traverse_to_depth(key, children, depth, cur_depth=1):
    """Traverse a list of children to a specified depth"""
    if cur_depth < depth:
        return concat(traverse_to_depth(key+k, ch, depth, cur_depth=cur_depth+1)
                      for k, ch in children.items())
    else:
        return [(key, concat(map(get_leaves, children.values())))]


def split_levels(code):
    """Split an ATC code into its 5 levels"""
    return [
        code[0],
        code[1:3],
        code[3],
        code[4],
        code[5:]
    ]

def get_level(code, level):
    """Get an ATC code up to the specified level"""
    return ''.join(split_levels(code)[:level])
