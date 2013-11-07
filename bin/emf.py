import dendropy
# EMF just for iterating over the file
# Push annotations to each node
__all__ = [ "EMF" ]

class EMF:
    def __init__(self, fn):
        self.fn = fn
        self._handle = open(fn)
        self._last = None

    def __iter__(self):
        return self

    def next(self):
        seqs = {}
        treedata = []

        if self._last:
            seqs.append(self._last)
            self._last = None
        
        for l in self._handle:
            if l == '':
                raise StopIteration

            if l == "\n":
                continue

            if not l.startswith("SEQ"):
                if not l == "DATA\n":
                    print l
                    raise ValueError("Malformed EMF file")
                else:
                    break
            f = l.rstrip().split()
            seqs[f[2]] = f


        for l in self._handle:
            if l.startswith("//"):
                break

            treedata.append(l)

        tree = dendropy.Tree.get_from_string(''.join(treedata), schema='newick')
        
        for l in tree.leaf_iter():
            f = seqs[str(l.taxon)]
            
            l.annotations.add_new('species', f[1])
            l.annotations.add_new('gene_id', f[7])

        return tree
        

# Test code
if __name__ == "__main__":
    emf = EMF("/Users/greg/Downloads/Compara.71.protein.nhx.emf")

    # it = iter(emf)
    # print type(it.next())

    for t in emf:
        print t.leaf_nodes()[0].annotations.find(name="species").value
