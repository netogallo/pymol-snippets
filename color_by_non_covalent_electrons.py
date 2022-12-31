import pymol.cmd as cmd
import pymol.selector as selector

def color_by_non_covalent_electrons(selection):
    """
    Color the carbon atoms that are covalently bound to an element
    which has electrons capable of non-covalent interactions.
    """
    
    selection = selector.process(selection)
    model = cmd.get_model(selection)

    color_1 = "color_1"
    color_2 = "color_2"
    color_3 = "color_3"
    cmd.set_color(color_1, [0.7,0,0])
    cmd.set_color(color_2, [0.8,0,0])
    cmd.set_color(color_3, [0.9,0,0])
    
    for (ix, atom) in [(ix, atom) for (ix, atom) in enumerate(model.atom) if atom.symbol == 'C']:
        bonds = [bond for bond in model.bond if ix in bond.__dict__['index']]

        for bond in bonds:
            [ix_bound] = [i for i in bond.__dict__['index'] if i != ix]
            atom_bound = model.atom[ix_bound]

            if atom_bound.symbol == 'O':
                cmd.color(color_2, "%s and index %s" % (selection, atom.index))
                # return "%s and index %s" % (selection, atom.index)
                break
                
