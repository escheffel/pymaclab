
# Author: Pearu Peterson
# Created: April 2011

__all__ = ['load_stoic_from_sbml',
           'load_stoic_from_text',
           ]

import re
import os
from collections import defaultdict

from .utils import obj2num


def load_stoic_from_sbml(file_name, split_bidirectional_fluxes=False):
    """ Return stoichiometry information of a network described in a SBML file.

    Parameters
    ----------
    file_name : str
      Path to SMBL file.

    split_bidirectional_fluxes : bool
      When True the bidirectional fluxes are split into two unidirectional fluxes.
      For example, the system ``A<=>B`` is treated as ``A=>B and B=>A``.

    Returns
    -------
    matrix : dict
      A stoichiometry matrix defined as mapping {(species, reaction): stoichiometry}.
    species : list
      A list of species names.
    reactions : list
      A list of reaction names.
    species_info : dict
    reactions_info : dict
    """
    from lxml import etree
    for ext in ['', '.xml', '.xml.gz', '.gz']:
        if os.path.isfile(file_name+ext):
            file_name = file_name+ext
            break
    tree = etree.parse(file_name)
    root = tree.getroot()
    assert root.tag.endswith ('sbml'), `root.tag`
    version = int(root.attrib['version'])
    level = int(root.attrib['level'])
    if level in [2,3]:
        default_stoichiometry = '1'
    else:
        default_stoichiometry = None
    compartments = {}
    species = []
    modifiers = []
    species_all = []
    reactions = []
    species_reactions = defaultdict (lambda:[])
    reactions_species = defaultdict (lambda:[])
    reactions_info = defaultdict(lambda:dict(modifiers=[],reactants=[],products=[],
                                            boundary_specie_stoichiometry={},annotation=[],
                                            compartments = set()))
    species_info = defaultdict(lambda:dict())
    matrix = {}
    for model in root:
        for item in model:
            if item.tag.endswith('listOfCompartments'):
                for compartment in item:
                    compartments[compartment.attrib['id']] = compartment.attrib
            elif item.tag.endswith('listOfSpecies'):
                for specie in item:
                    species_all.append(specie.attrib['id'])
                    species_info[specie.attrib['id']]['compartment'] = specie.attrib['compartment']
                    species_info[specie.attrib['id']]['name'] = specie.attrib.get('name', specie.attrib['id'])
            elif item.tag.endswith('listOfReactions'):
                for reaction in item:

                    reversible =eval(reaction.attrib.get('reversible', 'False').title())
                    reaction_id = reaction.attrib['id']
                    name = reaction.attrib.get('name', reaction_id)
                    assert reaction_id not in reactions,`reaction_id`
                    reactions.append(reaction_id)
                    reaction_index = len(reactions)-1
                    reactions_info[reaction_id]['name'] = name

                    if split_bidirectional_fluxes and reversible:
                        reaction_id2 = '%s_r' % (reaction_id)
                        assert reaction_id2 not in reactions,`reaction_id2`
                        reactions.append(reaction_id2)
                        reaction_index2 = len(reactions)-1
                        reactions_info[reaction_id2]['name'] = name+'_r'
                        reactions_info[reaction_id]['reversible'] = False
                        reactions_info[reaction_id2]['reversible'] = False
                    else:
                        reaction_id2 = reaction_index2 = None
                        reactions_info[reaction_id]['reversible'] = reversible

                    for part in reaction:
                        if part.tag.endswith ('listOfReactants'):
                            for reactant in part:
                                assert reactant.tag.endswith('speciesReference'), `reactant.tag`
                                specie_id = reactant.attrib['species']
                                stoichiometry = -obj2num(reactant.attrib.get('stoichiometry', default_stoichiometry))
                                reactions_info[reaction_id]['reactants'].append(specie_id)
                                try:
                                    specie_index = species.index(specie_id)
                                except ValueError:
                                    species.append(specie_id)
                                    specie_index = len(species)-1
                                assert stoichiometry,`stoichiometry`
                                matrix[specie_index, reaction_index] = stoichiometry
                                species_reactions[specie_index].append(reaction_index)
                                reactions_species[reaction_index].append(specie_index)                                
                                reactions_info[reaction_id]['compartments'].add(species_info[specie_id]['compartment'])
                                if reaction_index2 is not None:
                                    reactions_info[reaction_id2]['reactants'].append(specie_id)
                                    matrix[specie_index, reaction_index2] = -stoichiometry
                                    species_reactions[specie_index].append(reaction_index2)
                                    reactions_species[reaction_index2].append(specie_index)                                
                                    reactions_info[reaction_id2]['compartments'].add(species_info[specie_id]['compartment'])
                        elif part.tag.endswith ('listOfProducts'):
                            for product in part:
                                assert product.tag.endswith('speciesReference'), `product.tag`
                                specie_id = product.attrib['species']
                                stoichiometry = obj2num(product.attrib.get('stoichiometry', default_stoichiometry))
                                reactions_info[reaction_id]['products'].append(specie_id)
                                try:
                                    specie_index = species.index(specie_id)
                                except ValueError:
                                    species.append(specie_id)
                                    specie_index = len(species)-1
                                assert stoichiometry,`stoichiometry`
                                matrix[specie_index, reaction_index] = stoichiometry
                                species_reactions[specie_index].append(reaction_index)
                                reactions_species[reaction_index].append(specie_index)
                                reactions_info[reaction_id]['compartments'].add(species_info[specie_id]['compartment'])
                                if reaction_index2 is not None:
                                    reactions_info[reaction_id2]['products'].append(specie_id)
                                    matrix[specie_index, reaction_index2] = -stoichiometry
                                    species_reactions[specie_index].append(reaction_index2)
                                    reactions_species[reaction_index2].append(specie_index)
                                    reactions_info[reaction_id2]['compartments'].add(species_info[specie_id]['compartment'])
                        elif part.tag.endswith ('listOfModifiers'):
                            for modifier in part:
                                assert modifier.tag.endswith('modifierSpeciesReference'), `modifier.tag`
                                specie_id = product.attrib['species']
                                reactions_info[reaction_id]['modifiers'].append(specie_id)
                                reactions_info[reaction_id]['compartments'].add(species_info[specie_id]['compartment'])
                            continue
                        elif part.tag.endswith ('annotation'):
                            reactions_info[reaction_id]['annotation'].append(part.text)
                            continue
                        elif re.match(r'.*(kineticLaw|notes)\Z', part.tag):                            
                            continue
                        else:
                            print 'get_stoichiometry:warning:unprocessed reaction element: %r' % (part.tag)
                            continue


            elif re.match (r'.*(annotation|notes|listOfSpeciesTypes|listOfUnitDefinitions)\Z', item.tag):
                pass
            else:
                print 'get_stoichiometry:warning:unprocessed model element: %r' % (item.tag)

    return matrix, species, reactions, species_info, reactions_info
    

def load_stoic_from_text(text, split_bidirectional_fluxes=False):
    """ Parse stoichiometry matrix from a string.

    Parameters
    ----------
    text : str
      A multiline string where each line contains a chemical reaction
      description. The description must be given in the following
      form: ``<sum of reactants> (=> | <=) <sum of producats>``. For example,
      ``A + 2*B => C``. Lines starting with ``#`` are ignored.

    split_bidirectional_fluxes : bool
      When True the bidirectional fluxes are split into two unidirectional fluxes.
      For example, the system ``A<=>B`` is treated as ``A=>B and B=>A``.

    Returns
    -------
    matrix_data : dict
      A dictionary representing a stoichiometry matrix.

    species : list
      A list of row names.

    reactions : list
      A list of column names.

    species_info : dict
    reactions_info : dict
    """
    #TODO: fill up species_info and reactions_info dictionaries

    def _split_sum (line):
        for part in line.split('+'):
            part = part.strip()
            coeff = ''
            while part and part[0].isdigit():
                coeff += part[0]
                part = part[1:].lstrip()
            if not coeff:
                coeff = '1'
            if not part:
                continue
            yield part, eval (coeff)

    matrix = {}
    reactions = []
    species = []
    reactions_info = defaultdict(lambda:dict(modifiers=[],reactants=[],products=[],
                                            boundary_specie_stoichiometry={},annotation=[],
                                            compartments = set()))
    species_info = defaultdict(lambda:dict())
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith ('#'): continue
        if ':' in line:
            reaction_name, line = line.split (':',1)
            reaction_name = reaction_name.strip()
            line = line.strip()
        else:
            reaction_name = None
        if not line or line.startswith ('#'): continue

        reaction_string = line

        reversible = False
        left, right = line.split ('=')
        direction = '='
        if right.startswith('>'):
            right = right[1:].strip()
            direction = '>'
            if left.endswith ('<'):
                left = left[:-1].strip()
                reversible = True
        elif left.endswith ('>'):
            left = left[:-1].strip()
            direction = '>'
        elif left.endswith ('<'):
            left = left[:-1].strip()
            direction = '<'
        elif right.startswith ('<'):
            right = right[1:].strip()
            direction = '<'

        left_specie_coeff = list(_split_sum(left))
        right_specie_coeff = list(_split_sum(right))
        left_specie_names = [n for n,c in left_specie_coeff if n]
        right_specie_names = [n for n,c in right_specie_coeff if n]

        fname = ['R']
        rname = ['R']
        name0 = ''.join(left_specie_names)
        name1 = ''.join(right_specie_names)
        if name0:
            rname.append (name0)
        if name1:
            fname.append (name1)
            rname.append (name1)
        if name0:
            fname.append (name0)

        if direction=='<':
            if not reaction_name:
                reaction_name = '_'.join(fname)
                reaction_name2 = '_'.join(rname)
            else:
                reaction_name2 = 'r'+reaction_name
                if split_bidirectional_fluxes:
                    reaction_name = 'f'+reaction_name
        else:
            if not reaction_name:
                reaction_name2 = '_'.join(fname)
                reaction_name = '_'.join(rname)
            else:
                reaction_name2 = 'r'+reaction_name
                if split_bidirectional_fluxes:
                    reaction_name = 'f'+reaction_name

        reactions.append (reaction_name)
        reaction_index = reactions.index (reaction_name)
        if split_bidirectional_fluxes and reversible:
            reactions.append (reaction_name2)
            reaction_index2 = reactions.index (reaction_name2)
        else:
            reaction_index2 = None

        def matrix_add (i,j,c):
            v = matrix.get ((i,j))
            if v is None:
                matrix[i, j] = c
            else:
                matrix[i, j] = v + c

        for specie, coeff in left_specie_coeff:
            if specie not in species:
                species.append (specie)
            specie_index = species.index (specie)
            if direction=='<':
                if reaction_index2 is not None:
                    matrix_add(specie_index, reaction_index2, -coeff)
                matrix_add(specie_index, reaction_index, coeff)
            else:
                if reaction_index2 is not None:
                    matrix_add(specie_index, reaction_index2, coeff)
                matrix_add(specie_index, reaction_index, -coeff)
            
        for specie, coeff in right_specie_coeff:
            if specie not in species:
                species.append (specie)
            specie_index = species.index (specie)
            if direction=='<':
                if reaction_index2 is not None:
                    matrix_add(specie_index, reaction_index2, coeff)
                matrix_add(specie_index, reaction_index, -coeff)
            else:
                if reaction_index2 is not None:
                    matrix_add(specie_index, reaction_index2, -coeff)
                matrix_add(specie_index, reaction_index, coeff)

        if split_bidirectional_fluxes:
            reactions_info[reaction_name]['reversible'] = False
            reactions_info[reaction_name]['reactants'] = left_specie_names
            reactions_info[reaction_name]['products'] = right_specie_names
            reactions_info[reaction_name]['forward'] = reaction_name
            reactions_info[reaction_name]['reverse'] = None

            if reversible:
                reactions_info[reaction_name2]['reversible'] = False
                reactions_info[reaction_name2]['reactants'] = right_specie_names
                reactions_info[reaction_name2]['products'] = left_specie_names
                reactions_info[reaction_name2]['forward'] = reaction_name2
                reactions_info[reaction_name2]['reverse'] = None
            # TODO: set  reactions_info[reaction_name]['name']
        else:
            reactions_info[reaction_name]['reversible'] = reversible
            reactions_info[reaction_name]['reactants'] = left_specie_names
            reactions_info[reaction_name]['products'] = right_specie_names
            if reversible:
                reactions_info[reaction_name]['forward'] = 'f'+reaction_name
                reactions_info[reaction_name]['reverse'] = 'r'+reaction_name
            else:
                reactions_info[reaction_name]['forward'] = 'f'+reaction_name
                reactions_info[reaction_name]['reverse'] = None

            reactions_info[reaction_name]['name'] = reaction_string

    return matrix, species, reactions, species_info, reactions_info
