# -*- python -*-
#
#  This file is part of MS-DAS software
#
#  Copyright (c) 2014 - EBI-EMBL
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>, Marti Bernardo Faura
#  bernardo@ebi.ac.uk
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#
##############################################################################
from msdas.psites import PSites

__all__ = ["Requires", "SequenceTools", "Modification"]

# just a simple class to mimic behaviour of easydev.decorators.requires decorator
# this issue with the decorator is that it prevents sphinx to see the parametes/arguments

class Requires(object):
    """A fake decorator to require attribute to be present in a class instance.

    .. warning:: this is for developers only.

    This is not a decorator but a method to be called at the beginning of
    another method to check for the presence of an attribute.

    """
    def __init__(self):
        pass

    def _requires(self, attributes, msg):
        """Equivalent of easydev.decorators.requires

        This function can be used by a method to check that the class
        has the relevant attributes.

        """
        if isinstance(attributes, str):
            attributes = [attributes]
        elif isinstance(attributes, list):
            pass
        else:
            raise TypeError("attributes parameter must be a list of strings or a string")

        for attribute in attributes:
            if hasattr(self, attribute) == False:
                raise AttributeError("Missing attribute ({}) in the class. {}".format(attribute, msg))


class SequenceTools(object):
    """Some common tools for developers

    .. warning:: this is for developers only.

    """

    def _interpret_fasta_header(self, header):
        res = {}

        header = header.strip()
        _, entry_accession, entry_name = header.split("|")

        try:
            name, organism = entry_name.split("_")
        except:
            name, organism = (entry_name, "undefined")
        iso = entry_accession.split("-")
        assert len(iso) in [1,2]
        if len(iso)==2:
            iso = entry_accession.split("-")
            iso = iso[1]
            name_iso = name + "-" + iso
        else:
            name_iso = name
        res["accession"] = entry_accession
        res["name"] = name_iso
        res["organism"] = organism
        return res

    def get_protein_from_fasta_header(self, name):
        """Returns protein name given the full accession name

        ::

            >>> from msdas.tools import SequenceTools
            >>> t = SequenceTools()
            >>> name = 'sp|Q8IYB3-2|SRRM1_HUMAN'
            >>> res = t.get_protein_from_fasta_header(name)
            >>> res
            'SRRM1-2'

        """
        res = self._interpret_fasta_header(name)
        return res["name"]

    def get_proteinId_from_fasta_header(self, name):
        """Returns the accession name


        ::

            >>> from msdas.tools import SequenceTools
            >>> t = SequenceTools()
            >>> name = 'sp|Q8IYB3-2|SRRM1_HUMAN'
            >>> res = t.get_protein_from_fasta_header(name)
            >>> res
            'Q8IYB3-2'

        """
        res = self._interpret_fasta_header(name)
        return res["accession"]


class Modification(object):
    """Tools to manipulate the modification string found in mass spec data.


    Some Mass Spec data contains modifications encoded as follows::

        [3] S+79.9663|[5] S+79.9663

    For this sequence::

        RESPSPAPKPR

    Meaning there are 2 phospho sites  (separated by the pipe) at position 3 (an S)
    and 5 (another S). This class provides tools to handle the **modification**.

    ::

        m = Modification()
        m.get_individual_psites_from_modification("[3] S+79.9663|[5] S+79.9663")

    """
    _error_message = "Invalid modification. Use pipe character to separate them e.g., [6] S+79.663|[7] S+79.663"
    example = "[3] S+79.9663|[5] S+79.9663"

    def __init__(self, modif=None):
        """.. rubric:: construcor


        :param modif: a valid modification string. If not provided, can be
            provided later as an attribute

        """
        self._modification = modif
        self.psite_factory = PSites()
        # location are counted based on the peptide sequence rather than the full sequence
        self.mode = "relative" # TODO not used ??

    def _get_psites(self):
        modif = self.get_individual_psites()
        psites = self.psite_factory.sorted(modif)
        return psites
    psites = property(_get_psites, doc="")

    def _get_modification(self):
        return self._modification
    def _set_modification(self, modif):

        #if self.isvalid(modif) == False:
        #    raise ValueError(self._error_message)
        self._modification = modif
    modification = property(_get_modification, _set_modification)

    def _get_modif(self, modif=None):
        """

        If modif is None, return the attribute, which is possibly set to
        a valid modification or set to None (default).

        If modif provided is not None, it is checked for validity.


        """
        if modif == None:
            if self.modification == None:
                return
            else:
                return self.modification
        else:
            #if self.isvalid(modif) == True:
            return modif
            #else:
            #    raise ValueError(self._error_message)

    def get_individual_psites(self, modif=None):
        """Returns individual psite location given a full string of modifications

        ::

            >>> from msdas import Modification
            >>> m = Modification()
            >>> m.get_individual_psites("[5] S+7777| [1] T+5555")
            ['S5', 'T1']

        """
        modif = self._get_modif(modif)

        # split according to pipe character if any
        modif = modif.split("|")

        # let us get rid of the + followed by numbers
        modif = [x.split("+")[0] for x in modif]

        # let us get rid of spaces, and brackets
        modif = [x.replace(" ","") for x in modif]
        modif = [x.replace("[","") for x in modif]
        modif = [x.replace("]","") for x in modif]

        # psites should be letter + location
        # Each modif should be correct that is a digit ending with a valid
        #psite S,Y,T
        modif = [x[-1] + x[0:-1] for x in modif]

        return modif

        # order

    def psites2modif(self, psites):
        """Transform a psites to a modification-like string

        :param str psites: e.g. S1^T2 see :class:`msdas.psites.PSites` for details.

        Not implemented if + characters are present in the psites string

        """
        self.psite_factory.isvalid(psites)

        AND = self.psite_factory.AND
        psites = psites.split(AND)
        psites = " | ".join(["["+x[1:]+"]"+x[0] for x in psites])
        return psites

    def modif2psites(self, modif=None):
        """Return single psites stirng from modification string


        ::

             >>> from msdas import Modification
             >>> m = Modification()
             >>> m.modif2psites("[5] S+7777| [1] T+5555")
             'T1^S5'

        .. note:: psites are sorted by position

        .. seealso:: :meth:`get_individual_psites`
        """
        modif = self._get_modif(modif)
        modif = self.get_individual_psites(modif)
        psites = self.psite_factory.sorted(modif)
        return psites






