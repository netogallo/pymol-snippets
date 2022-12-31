import pymol.cmd as cmd
import pymol.selector as selector

def color_by_non_covalent_electrons(selection):
    """
    Color the carbon atoms that are covalently bound to an element
    which has lots of electrons that can interact non-covalently. The
    stronger the red color, the more electrons are present in the
    site.
    """
    
    selection = selector.process(selection)
    cmd.color('green', selection)
    model = cmd.get_model(selection)

    color_3 = "color_1"
    color_4 = "color_2"
    color_5 = "color_3"
    color_6 = "color_4"
    cmd.set_color(color_3, [0.9,0.9,0])
    cmd.set_color(color_4, [0.9,0.6,0])
    cmd.set_color(color_5, [0.9,0.3,0])
    cmd.set_color(color_6, [0.9,0,0])

    carbons = [(ix, atom) for (ix, atom) in enumerate(model.atom) if atom.symbol == 'C']
    carbons_ix = set([i for (i,_) in carbons])

    bonds_dict = {}

    electrons_outter_shell = { 'O' : 6, 'N' : 5 }

    for bond in model.bond:
        for ix in bond.__dict__['index']:
            if not ix in bonds_dict:
                bonds_dict[ix] = [bond]
            else:
                bonds_dict[ix].append(bond)
    
    def search_non_covalent(start_ix):
        
        visited = set(carbons_ix)
        queue = [start_ix]
        non_covalent_count = 0

        while(len(queue) > 0):
            ix = queue.pop()
            visited.add(ix)
            atom = model.atom[ix]

            for bond in bonds_dict[ix]:
                for i in bond.__dict__['index']:
                    if not i in visited:
                        queue.append(i)

            electrons = electrons_outter_shell.get(atom.symbol)

            if electrons:
                total_bonds = sum([bond.__dict__['order'] for bond in bonds_dict[ix]])
                non_covalent_count = non_covalent_count + electrons - total_bonds

        return non_covalent_count
     
    for (ix, atom) in [(ix, atom) for (ix, atom) in enumerate(model.atom) if atom.symbol == 'C']:
        bonds = bonds_dict[ix]
        score = 0

        for bond in bonds:
            ix_bounds = [i for i in bond.__dict__['index'] if i not in carbons_ix]
            score = score + (search_non_covalent(ix_bounds[0]) if len(ix_bounds) > 0 else 0)

        if score > 6:
            cmd.color(color_6, "%s and index %s" % (selection, atom.index))
        if score > 5:
            cmd.color(color_5, "%s and index %s" % (selection, atom.index))
        elif score > 4:
            cmd.color(color_4, "%s and index %s" % (selection, atom.index))
            # return "%s and index %s" % (selection, atom.index)
        elif score > 3:
            cmd.color(color_3, "%s and index %s" % (selection, atom.index))
#        elif score > 2:
#            cmd.color(color_2, "%s and index %s" % (selection, atom.index))                

# color_by_non_covalent_electrons('6lgg')

pymol.cmd.extend('color_by_non_covalent_electrons', color_by_non_covalent_electrons)
