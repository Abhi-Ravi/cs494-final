import pubchempy as pubchem


def checkLipinski(compound_input):
    compounds = pubchem.get_compounds(compound_input, 'name')

    c = pubchem.Compound.from_cid(compounds[0].cid)

    if c.h_bond_donor_count > 5:
        return False

    if c.h_bond_acceptor_count > 10:
        return False

    if c.exact_mass > 500:
        return False

    if c.xlogp > 5:
        return False

    return True


compound_input = input()

if checkLipinski(compound_input):
    print("Passes Lipinski rule of five")
else:
    print("Fails Lipinski rule of five")