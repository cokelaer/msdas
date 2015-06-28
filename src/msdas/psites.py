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
import re
from easydev import Logging

__all__ = ['PSites']


class PSites(Logging):
    """Utilities to manipulate and validate phosphorylation sites (psites) encoding

    Psites stands for phosphorylation sites. There are encoded in different w
    ays depending on the input data files. This class
    provides methods to help decipher the different format and convert them
    into a common one.

    Psites indicate the position of a phosphorylation on a given sequence of
    amino acids as well as the amino acid on which the phosphorylation occured. For
    instance::

        S20

    means there is a phosphorylation on Serine at poistion 20 of the protein sequence.
    Possibly, you may have more than one phosphorylation. In order to encode several
    locations, we use the **^** character to represent an **AND**. Mote that
    sometimes a location may be ambiguous. In such case, we use the **+**  character
    to represent an **OR**.

    We use the following conventions.

    #. no spaces
    #. aminod acids are encoded with capitalised letters
    #. S25^T40 means there are 2 psites at position 25 AND 40
    #. S25+S26 means there is an ambiguity and the psite is located at poisiton 25 OR 40
    #. S25+S26^T40 means there are 2 psites, one of which is ambiguous (S25).

    The number of psites is equal to the number of ^ characters plus one.

    """
    OR = "+"
    AND = "^"

    def __init__(self, verbose=True):
        super(PSites, self).__init__(level=verbose)
        # valid phosphos are S, T, Y; M is an oxidation that may not be required
        self.valid_letters = ['S', 'T', 'Y', 'M']


    def append_oxidation(self, sequence, psites):
        """Update psite string with Oxydation location if present in the sequence

        Sometimes, the tag Oxydation is present in the sequence but not reflected
        in the psite string. This function update the psite according to the
        Oxidation tags.

        :param str sequence: sequence of amino acid with position of phospho and oxidation
            as in the following example: KVS(Phospho)VGSM(Oxidation)G
        :param str sites: a valid psite string (capitalised letters)
        :return: Updated psite string with oxidation locations. See example below.


        .. doctest::

            >>> from msdas.psites import PSites
            >>> p = PSites()
            >>> sequence = "KVS(Phospho)VGSM(Oxidation)GS(Phospho)GK"
            >>> psites = "S956^S962"
            >>> p.append_oxidation(sequence, psites)
            'S956^M960^S962'

        .. todo:: somehow we do not need the psite string as a parameter. Can be extracted
            from the sequence anyway.

        """
        if self.isvalid(psites) == False:
            raise ValueError("invalid sites provided () must be separated by ^ character" % psites)

        if "Oxidation" not in sequence:
            return psites
        psites = psites.split(self.AND)
        absolute_position = int(psites[0][1:]) - sequence.replace("(Oxidation)","").index("(Phospho)")

        for occ in re.finditer("(Oxidation)", sequence):
            seq = sequence[0:occ.start()].replace("(Phospho)", "").replace("(Oxidation)", "")
            seq = seq.replace("(", "")
            nucleo = seq[-1]
            new_psite =  nucleo + str(absolute_position + len(seq))
            psites.append(new_psite)

        # now, we sorte the psite
        psites = self.sorted(psites)
        return psites

    def get_unique_psites(self, psites_list):
        """

        handles both AND and OR characters
        """
        psites_list = [x for y in psites_list for x in y.split("^")]
        psites_list = [x for y in psites_list for x in y.split("+")]

        return list(set(psites_list))

    def get_factorised_psites(self, psites_list):
        """Given a list of psites, find common and ambiguous parts and return new string

        Only the AND psites are factorised.

            >>> p = PSites()
            >>> p.get_factorised_psites(['S177^Y182', 'T180^Y182'])
            'S177+T180^Y182'
            >>> p.get_factorised_psites(["S272^T277^S279", "S272^S275^S279"])
            'S275+T277^S272^S279'

        The returned string in the last example tells us that there are 3
        phosphosites (2 ^ character + 1). The first psiste is ambiguous (S275+T277)

        Note this example::

            >>> p.get_factorised_psites(["S1+S4^S3", "S2+S1^S3"])
            'S1+S1+S2+S4^S3'

        .. warning:: in this example, be aware that (i), even though psites are
            not ordered by location in the input, there are in the output and
            (ii) if there are duplicates sites (e.g., S1)  there are not
            simplified for now.

        """
        all_sites = set([y for x in psites_list for y in x.split("^")])

        # this is the set of common psites across all psites
        common = self.get_common_psites(psites_list)

        # and the ones that are not
        not_common = all_sites.difference(common)
        # sorted concatenate with ^; let us replace by +
        not_common = self.sorted(list(not_common)).replace("^", "+")

        psites = [not_common] + list(common)
        psites = [psite for psite in psites if psite] # remove empty strings
        if len(psites)>1:
            psites= "^".join([not_common]+list(common))

        psites = self.sorted(psites.split("^"))

        return psites

    def get_common_psites(self, psites_list):
        """Given a list of psites of the form Y182^T180, returns common factor


            >>> p = PSites()
            >>> p.get_common_psites(['S177^Y182', 'T180^Y182'])
            {'Y182'}

        The OR characters are ignored.

        .. note:: used by :meth:`get_factorised_psites`

        """
        all_sites = set([y for x in psites_list for y in x.split(self.AND)])
        sets = [set(psites.split(self.AND)) for psites in psites_list]
        # note the * in front of the list of sets to make it a list of arguments
        return all_sites.intersection(*sets)

    def isvalid(self, psites):
        """Return true is site is using a valid format

        Checks that:

        #. letters are capitalised and in the :attr:`valid_letters` (default T, Y, S, M)
        #. there is no spaces
        #. number of aminod acid is compatible with separators (+/^). e.g, T4T5 is not
           correct. It should be either S4+S5 or S4^S5.

        :param str psite: a psites string (e.g., S2^T5)
        :return: boolean

        .. seealso:: the convention is explained in the class documentation.

        ::

            >>> from msdas import PSites
            >>> p = PSites()
            >>> p.isvalid("T1+S2^T")
            False
        """
        if " " in psites:
            self.logging.warning("Invalid psite with spaces inside")
            return False

        psites = [x for y in psites.split(self.OR) for x in y.split(self.AND)]
        # check letters
        for psite in psites:
            found = list(re.finditer("[%s]" % ",".join(self.valid_letters), psite))
            if len(found)==0:
                self.logging.warning("Found an invalid letter (must be captialised and in {}). You provided {},{}".format(self.valid_letters, psite, psites))
                return False
        # check that there are number after letter
        for psite in psites:
            if len(psite)==1:
                self.logging.warning("Invalid psite missing position in {}, {}?".format(psite,psites))
                return False
            number = psite[1:]
            if number.isdigit() == False:
                self.logging.warning("Invalid location; should be numeric after Phospho {},{}?".format(psite,psites))
                return False

        return True

    def remove_duplicated(self, psites):
        """Remove duplicates psites in a psite string


         .. doctest::

            >>> from msdas.psites import PSites
            >>> p = PSites()
            >>> p.remove_duplicated("S3^S4^S3")
            'S3^S4'

        """
        psites = psites.split(self.AND)
        psites = list(set(psites))
        psites = self.sorted(psites)
        return psites

    def remove_spaces(self, psite):
        """Removes spaces found in a psite string

        Sometimes, psites are separated by ;. These characters are kept
        in this function but should be replaced once the meaning is known
        by AND or OR character.

        """
        psite = psite.strip()
        psite = ";".join([x.strip() for x in psite.split(";")])
        return psite

    def sort_psites_ors_only(self, psites):
        """Sort psites

        :param str psites: a valid psites string.
        :return: modified psites

        ::

            >>> from msdas import PSites
            >>> p = PSites()
            >>> p.sort_psites_ors_only("S2+S1")
            'S1+S2'

        Psites that contain AND character are not processed.

        """
        if self.AND in psites:
            return psites
        else:
            f = lambda x,y: cmp(int(x[1:].split(self.OR)[0]), int(y[1:].split(self.OR)[0]))

            return "+".join(sorted(psites.split("+"), cmp=f))


    def sorted(self, psites_list):
        """Sort the psites found in the psites list using their positions

        If there is an ambiguous site, the + is ignored and only the first psite
        is used for the position.

        .. doctest::

            >>> from msdas.psites import PSites
            >>> p = PSites()
            >>> p.sorted(["S3", "S1", "S4+S5"])
            'S1^S3^S4+S5'
            >>> p.sorted(["S3", "S1", "S8+S7"])
            'S1^S3^S7+S8'

        """

        # check validity
        f = lambda x,y: cmp(int(x[1:].split(self.OR)[0]), int(y[1:].split(self.OR)[0]))

        # First, we need to order internally the psites that contain ORs
        # e.g., S8+S7 must be re-ordered as S7+S8
        for psite in psites_list:
            if self.isvalid(psite) == False:
                self.logging.error("Invalid psite %s provided" % psite)

        # each item that has an OR will be split, sorted and joined back
        psites_list = [self.OR.join(sorted(psite.split(self.OR), cmp=f)) for psite in psites_list]

        # now, we do the actual sorting of the list
        psites =  sorted(psites_list, cmp=f)
        psites = self.AND.join(psites)
        return psites

    def remove_duplicated(self, psites):
        """Remove duplicates psites in a psite string


         .. doctest::

            >>> from msdas.psites import PSites
            >>> p = PSites()
            >>> p.remove_duplicated("S3^S4^S3")
            'S3^S4'


        """
        psites = psites.split(self.AND)
        psites = list(set(psites))
        psites = self.sorted(psites)
        return psites
