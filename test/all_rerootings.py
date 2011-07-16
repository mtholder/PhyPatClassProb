#!/bin/env python
import sys
import dendropy
stream = open(sys.argv[1], 'rU')
tree = dendropy.TreeList(stream=stream, schema="nexus")[0]
node_list = [nd for nd in tree.preorder_node_iter()]
node_list.pop(0)
print(tree.as_string('newick'))
for nd in node_list:
    if not nd.is_leaf():
        tree.reroot_at_node(nd, update_splits=False)
        print(tree.as_string('newick'))
