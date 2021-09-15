def atom_radius(atom, library, fallback=1.50):
    """
    Function to return the atomic radius from the supplied library

    Input:
    atom = Atom symbol (up to 2 characters)
    library = Dictionary of atomic radii, with atom symbols as keys and radii as values
    fallback = if an atom is not in the dictionary, use this value as the fallback radius

    Output: returns either the radius from the library or the fallback value
    """
    if atom in library:
        return library[atom]
    else:
        return fallback
