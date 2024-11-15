# Copyright (C) 2018,2019 Luiz Felippe S. Rodrigues, Luke Chamandy
#
# This file is part of Magnetizer.
#
# Magnetizer is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Magnetizer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Magnetizer.  If not, see <http://www.gnu.org/licenses/>.
#
import numpy as np
import astropy.units as u


class BinningObject(object):
    """
    The attributes contain the mass bins and filters which can be applied to
    select specific mass bins.
    """
    def __init__(self, magnetizer_run, z=0.0, bins=None,
                 bin_array=None, extra_filter=None):

        self.run = magnetizer_run # Reference to run object (to avoid mistakes)
        if bins is None:
            if bin_array is None:
                raise ValueError
            M_bins = bin_array
            bins = []
            for m_m, m_p in zip(M_bins[:-1],M_bins[1:]):
                bins.append((m_m, m_p))

        self.redshift = z
        self.bins = tuple(bins)
        self.nbins = len(self.bins)
        self.masks = [None]*self.nbins
        self.bin_counts = [0]*self.nbins
        self._quantitytype = None
        self.extra_filter = extra_filter

    def _compute_bin_filter(self,quantity, bin_interval):
        select  = quantity >  bin_interval[0]
        select *= quantity <= bin_interval[1]
        return select

    def _update_masks(self,quantity):
        for i, interval in enumerate(self.bins):
            if self.masks[i] is None:
               self.masks[i] = self._compute_bin_filter(quantity, interval)
               if self.extra_filter is not None:
                  self.masks[i] *= self.extra_filter
            else:
               self.masks[i] *= self._compute_bin_filter(quantity, interval)
            self.bin_counts[i] = self.masks[i].sum()

    def __repr__(self):
        return_str = '[BinningObject]\n'
        return return_str + self.__str__()

    def __str__(self):
        return_str = '{} bins:\n'.format(self._quantitytype)
        for (b1, b2), bc in zip(self.bins, self.bin_counts):
            return_str += '[{0}, {1}] {2} galaxies\n'.format(b1, b2, bc)
        return return_str


class MassBinningObject(BinningObject):
    def __init__(self, magnetizer_run, z=0.0, bins=None, stellar_mass=True,
                 gas_mass=False, include_bulge=True, extra_filter=None,
                 bin_array=10**np.array([8.,8.75,9.5,10.25,11.])*u.Msun):

        super(MassBinningObject, self).__init__(magnetizer_run, z=z, bins=bins,
                               bin_array=bin_array, extra_filter=extra_filter)

        self.masks = [None]*len(self.bins)

        mass = None

        if stellar_mass:
            mass = magnetizer_run.get('Mstars_disk', z).copy()
            str_type = 'stellar'

            if include_bulge:
                mass += magnetizer_run.get('Mstars_bulge', z)
                str_disk = ' '
            else:
                str_disk = ' disc '

        if gas_mass:
            gas_mass = magnetizer_run.get('Mgas_disk', z).copy()
            if mass is None:
                mass = gas_mass
            else:
                mass += gas_mass
            str_type = 'gas'

        if stellar_mass and gas_mass:
            str_type = 'total'

        self._quantitytype = str_type + str_disk + 'mass'
        self._update_masks(mass)



