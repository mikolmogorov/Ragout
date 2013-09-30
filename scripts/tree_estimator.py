from collections import namedtuple
import math

#TreeNode = namedtuple("TreeNode", ["value", "left", "right"])
Edge = namedtuple("Edge", ["length", "node"])

class TreeNode:
    def __init__(self, values=[], left=None, right=None):
        self.values = values
        self.left = left
        self.right = right


def build_tree():
    LEN = 0.01
    leaves = [TreeNode(["a"]) for a in xrange(4)]
    leaves[0].values = ["b"]
    leaves[1].values = ["b"]
    #leaves[3].values = ["c"]
    print map(lambda l : l.values, leaves)
    v1 = TreeNode([], Edge(LEN, leaves[0]), Edge(LEN, leaves[1]))
    v2 = TreeNode([], Edge(LEN, leaves[2]), Edge(LEN, leaves[3]))
    root = TreeNode([], Edge(LEN, v1), Edge(LEN, v2))
    return root


def print_tree(root):
    def rec_print(root):
        print "{",
        if root.left:
            rec_print(root.left.node)
        print str(root.values),
        if root.right:
            rec_print(root.right.node)
        print "}",
    rec_print(root)
    print ""



def go_up(root):
    if not root.left and not root.right:
        return root.values

    left = [] if not root.left else go_up(root.left.node)
    right = [] if not root.right else go_up(root.right.node)
    intersect = [l for l in left if l in right]

    if intersect:
        root.values = intersect
    else:
        root.values = left + right
    return root.values


def go_down(root):
    if root.left:
        intersect = [l for l in root.left.node.values if l in root.values]
        if intersect:
           root.left.node.values = intersect
        go_down(root.left.node)

    if root.right:
        intersect = [l for l in root.right.node.values if l in root.values]
        if intersect:
           root.right.node.values = intersect
        go_down(root.right.node)


def subs_cost(v1, v2, length):
    MU = 0.01
    if v1 == v2:
        return -MU * length
    else:
        #print math.log(1 - math.exp(-MU * length))
        return math.log(1 - math.exp(-MU * length))


def tree_likelihood(root):
    if not root.left:
        left_prob = 0
    else:
        left_prob = (tree_likelihood(root.left.node) +
                    subs_cost(root.values[0], root.left.node.values[0], root.left.length))

    if not root.right:
        right_prob = 0
    else:
        right_prob = (tree_likelihood(root.right.node) +
                    subs_cost(root.values[0], root.right.node.values[0], root.right.length))
    #print left_prob, right_prob
    return left_prob + right_prob


def enumerate_trees(root):
    left = enumerate_trees(root.left.node) if root.left else [None]
    right = enumerate_trees(root.right.node) if root.right else [None]

    trees = []
    for val in root.values:
        for l in left:
            for r in right:
                el = Edge(root.left.length, l) if l else None
                er = Edge(root.right.length, r) if r else None
                trees.append(TreeNode([val], el, er))
    return trees


def max_likelihood_score(tree):
    go_up(tree)
    go_down(tree)
    print_tree(tree)
    max_score = float("-inf")
    for t in enumerate_trees(tree):
        likelihood = tree_likelihood(t)
        print_tree(t)
        print likelihood

        max_score = max(max_score, likelihood)
    return max_score


def main():
    tree = build_tree()
    #print_tree(tree)
    print max_likelihood_score(tree)


if __name__ == "__main__":
    main()
