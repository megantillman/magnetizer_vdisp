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
"""
Contains routines that compute auxiliary quantities.
This shall be updated later to include more complex properties.
"""
from numpy import pi, arctan, sqrt
import numpy as np
import re
import astropy.units as u
import astropy.constants as const

def compute_extra_quantity(qname, mag_run, gal_id=None, z=None, ivol=0,
                           cache=True):
    """
    Computes extra Magnetizer quantites.

    This function needs to be called with either gal_id or z different
    from None - i.e. the extra quantity qname will be computed for a given
    galaxy at various redshifts or for various galaxies at a single redshift.

    Some of the available quantities are:
    * '|Bp|'/'|Br|'/'|Bz|' - Absolute value of a B component
    * 'Btot' - Total (large scale) magnetic field strength
    * 'Bfull' - Total magnetic field strength
    * r'Bmax' - Maximum field strength
    * r'b' - Turbulent magnetic field strength
    * r'growth' - Growth rate of the B (assuming no radial variation)
    * r'growth_max' - Growth rate of the maximum field strength
    * r'Bfloor' - Target strength of the floor (independently estimated)
    * 'V' - Rotation curve
    * 'D' - Dynamo number
    * 'Dc'/'D_'Dcrit' - critical dynamo number
    * r'D_Dc' - D/Dc
    * 'R_u' - Reynolds number associated with vertical outflow
    * 'p' - Magnetic pitch angle
    * r'B_Beq' - |B|/B_eq
    * r'\tau\Omega'
    * 'q'
    * 'l/h'
    * 'h/r'

    Also, for convenience: any quantity, Q, can be computed at the radius of
    the maximum field strength [i.e. Q(rmax) where B(rmax)=Bmax]
    using 'Q_at_Bmax'.

    Parameters
    ----------
    qname : string
        code-name of the quantity

    mag_run : MagnetizerRun

    gal_id : integer
        index of the required galaxy. Having this different from None
        signals that the redshift evolution of a given galaxy should be
        computed
    z : float
        target redshift . Having this different from None signals that
        the quantity needs to be computed at this redshift only for a
        range of galaxies.
    cache : bool
        Caches any intermediate quantities requested for the calculation

    Returns
    -------
    astropy quantity
        if gal_id!=None, this will have shape (ngrid,nz) or (nz)
        if z!=None, this will have shape (ngals,ngrid) or (ngals)
    """

    # Makes sure the correct method is used
    if (gal_id is not None) and (z is None) :
        get = lambda quantity: mag_run.get_galaxy(quantity, gal_id, ivol=ivol)
    elif (gal_id is None) and (z is not None):
        get = lambda quantity, redshift=z: mag_run.get(quantity, redshift,
                                                       cache=cache)
    else:
        raise ValueError('Must choose either gal_id or z')

    if qname == 'V':
        quantity = get('Omega')*get('r')

    elif qname == 'Dgen':

        S = get('Shear')
        eta_t= get('etat')
        alpha = get('alp')
        h = get('h')

        quantity = alpha * S * h**3 / (eta_t**2)
        quantity = quantity.cgs

    elif qname == 'D':
        S = get('Shear')
        eta_t= get('etat')
        alpha = get('alp_k')
        h = get('h')

        quantity = alpha * S * h**3 / (eta_t**2)
        quantity = quantity.cgs

    elif qname == 'D_diag1':
        S = get('Shear')
        Om = get('Omega')
        h = get('h')
        try:
            v = get('v')
        except ValueError:
            print('Warning: v not in the output, using P_ISM_SOUND_SPEED_KM_S instead.')
            v = mag_run.parameters.ISM_and_disk['P_ISM_SOUND_SPEED_KM_S']
            v *= u.km/u.s

        quantity = 9. * h**2  * Om * S/ v**2
        quantity = quantity.cgs

    elif qname == 'D_diag2':
        S = get('Shear')
        h = get('h')
        l = get('l')

        try:
            v = get('v')
        except ValueError:
            print('Warning: v not in the output, using P_ISM_SOUND_SPEED_KM_S instead.')
            v = mag_run.parameters.ISM_and_disk['P_ISM_SOUND_SPEED_KM_S']
            v *= u.km/u.s

        quantity = 9. * h**3 * S /v /l**2
        quantity = quantity.cgs

    elif qname == 'D_diag3':
        S = get('Shear')
        h = get('h')

        try:
            v = get('v')
        except ValueError:
            print('Warning: v not in the output, using P_ISM_SOUND_SPEED_KM_S instead.')
            v = mag_run.parameters.ISM_and_disk['P_ISM_SOUND_SPEED_KM_S']
            v *= u.km/u.s

        quantity = 9. * h * S /v
        quantity = quantity.cgs

    elif qname == 'q_at_Bmax':
        quantity = -get('Shear_at_Bmax')/get('Omega_at_Bmax')

    elif qname == 'Bmax_b_at_Bmax':
        quantity = get('Bmax_Beq') / mag_run.parameters.dynamo['FMAG']

    elif qname == 'D1_D_at_Bmax':
        D1 = get('D_diag1_at_Bmax')
        D = get('D_at_Bmax')

        quantity = D1/D

    elif qname == 'D2_D_at_Bmax':
        D2 = get('D_diag2_at_Bmax')
        D = get('D_at_Bmax')

        quantity = D2/D

    elif qname == 'D3_D_at_Bmax':
        D3 = get('D_diag3_at_Bmax')
        D = get('D_at_Bmax')

        quantity = D3/D

    elif qname in (r'D_{{\rm crit}}','Dc','Dcrit'):
        Ru = get('R_u')
        Cu = 0.25 # hard-coded at input_constants module
        quantity = - (pi/2.)**5 * (1. + 4./pi**2 *Cu * Ru)**2
        quantity = quantity.cgs


    elif qname == 'alp_k':
        # Uses Krause formula to compute alpha_k if it is not in the output

        C_alp = mag_run.parameters.dynamo['C_ALP']
        l = get('l')
        h = get('h')
        try:
            v = get('v')
        except ValueError:
            print('Warning: v not in the output, using P_ISM_SOUND_SPEED_KM_S instead.')
            v = l/l
            v *= mag_run.parameters.ISM_and_disk['P_ISM_SOUND_SPEED_KM_S']
            v *= u.km/u.s
        Om = get('Omega')

        quantity = C_alp*l**2/h*Om
        quantity = quantity.to(u.km/u.s)
        quantity = np.minimum(quantity,v)


    elif qname == 'R_u':
        quantity = get('Uz') * get('h') / get('etat')

    elif qname == '|Bp|':
        quantity = np.abs(get('Bp'))

    elif qname == '|Br|':
        quantity = np.abs(get('Br'))

    elif qname == '|Bz|':
        quantity = get('Bzmod')

    elif qname == 'p':
        quantity = arctan(get('Br')/get('Bp'))
        quantity = quantity.to(u.degree)

    elif qname == 'q':
        quantity = -get('Shear')/get('Omega')

    elif qname == 'l/h':
        quantity = get('l')/get('h')

    elif qname == 'h/r':
        quantity = get('h')/get('r')
        quantity = quantity.cgs

    elif qname == r'\tau\Omega':
        quantity = get('tau')*get('Omega')

    elif qname == r'Btot':
        quantity = sqrt(get('Bp')**2 +
                    get('Br')**2 +
                    get('Bzmod')**2)

    elif qname == r'Bfull':
        quantity = sqrt(get('Btot')**2 + get('b')**2)

    elif qname == r'B_Beq':
        quantity = sqrt(get('Bp')**2 +
                    get('Br')**2 +
                    get('Bzmod')**2)
        quantity /= get('Beq')

    elif qname == r'Dgen_Dc':
        quantity = get('Dgen')/get('Dc')

    elif qname == r'D_Dc':
        quantity = get('D')/get('Dc')

    elif qname == r'Beq':
        n = get('n')*const.m_p
        v  = mag_run.parameters.ISM_and_disk['P_ISM_SOUND_SPEED_KM_S']
        v *= mag_run.parameters.ISM_and_disk['P_ISM_KAPPA']
        v *= u.km/u.s

        cgs_gauss = u.cm**-0.5 * u.g**0.5 /u.s

        quantity = np.sqrt(4*pi*n) * v
        quantity = quantity.cgs / cgs_gauss * 1e6*u.microgauss

    elif qname == r'Bfloor':
        r = get('r')
        h = get('h')
        Beq = get('Beq')
        #l = get('l') # This is not what is actually used in the code!
        l = mag_run.parameters.ISM_and_disk['P_ISM_TURBULENT_LENGTH']*u.pc
        kappa = mag_run.parameters.dynamo['P_FLOOR_KAPPA']

        Delta_r = l * kappa
        fmag = 0.5
        Ncells= np.abs(3.*r*Delta_r*h/l**3)
        brms= fmag*Beq

        quantity = np.exp(-Delta_r/2./r)*brms/Ncells**(0.5)*l/Delta_r/3.
        quantity = quantity.to(u.microgauss)

    elif qname == r'growth':
        if z is not None:
            # Redshift selection
            iz = mag_run._closest_redshift_idx(z)
            if iz==0:
                return 0.0 * get('r').base / u.Gyr

            prev_z = mag_run.redshifts[iz-1]

            Btot1 = sqrt(get('Bp')**2 +
                         get('Br')**2 +
                         get('Bzmod')**2)
            Btot0 = sqrt(get('Bp',redshift=prev_z)**2 +
                         get('Br',redshift=prev_z)**2 +
                         get('Bzmod',redshift=prev_z)**2)

            delta_lnB = np.log(Btot1/Btot0)
            delta_t = mag_run.times[iz] - mag_run.times[iz-1]

            quantity = delta_lnB/delta_t
        else:
            # Profile selection
            Btot = sqrt(get('Bp')**2 +
                         get('Br')**2 +
                         get('Bzmod')**2)

            quantity = np.empty_like(Btot.base)
            quantity[:,0] = 0.0*u.Gyr

            delta_lnB = np.log(Btot[:,1:]/Btot[:,:-1])
            delta_t = mag_run.times[1:] - mag_run.times[:-1]

            for prof in delta_lnB:
                prof /=delta_t

            quantity[:,1:] = delta_lnB

    elif qname == r'growth_max':
        if z is not None:
            # Redshift selection
            iz = mag_run._closest_redshift_idx(z)

            if iz==0:
                return np.nan * get('Mstars_disk').base / u.Gyr

            prev_z = mag_run.redshifts[iz-1]

            Bmax1 = get('Bmax')
            Bmax0 = get('Bmax',redshift=prev_z)

            delta_lnB = np.log(Bmax1/Bmax0)
            delta_t = mag_run.times[iz] - mag_run.times[iz-1]

            quantity = delta_lnB/delta_t
        else:
            # Profile selection
            Bmax = get('Bmax')

            quantity = np.empty_like(Bmax.base)/u.Gyr

            quantity[0] = 0.0

            delta_lnB = np.log(Bmax[1:]/Bmax[:-1])
            delta_t = mag_run.times[1:] - mag_run.times[:-1]

            quantity[1:] = delta_lnB/delta_t

    elif qname == r'Bmax':
        quantity = get('Btot')
        quantity = __get_profile_max(quantity, gal_id=gal_id, z=z)

    elif qname == r'rmax_rdisk':
        quantity = get('rmax')/get('r_disk')

    elif qname == r'hmax_rdisk':
        quantity = get('h_at_Bmax')/get('r_disk')
        quantity = quantity.cgs

    elif qname == r'hmax_rmax':
        quantity = get('h_at_Bmax')/get('rmax')
        quantity = quantity.cgs

    elif qname in ('BoT', 'B/T'):
        Mbulge = get('Mstars_bulge')
        quantity = Mbulge/(Mbulge + get('Mstars_disk'))
    elif qname == r'b':
        quantity = get('Beq').copy()
        quantity *= mag_run.parameters.dynamo['FMAG']

    elif qname == r'Bavg':
        if z is not None:
            Btot = get('Btot')
            r = get('r')
            # Area weigthed average:
            # Bavg = sum( B*2*pi*r*dr)/sum(2*pi*r*dr)
            # Cancelling constant terms
            # Bavg = sum(B*r)/sum(r)
            quantity = (Btot*r).sum(axis=1)/r.sum(axis=1) #Sums for each galaxy
        else:
            #TODO
            raise NotImplementedError

    elif qname == r'Beavg':
        if z is not None:
            Btot = get('Btot')
            r = get('r')
            h = get('h')
            # Volume weigthed average (B associated with average energy density):
            # Beavg = sqrt( 8pi*sum( B**2/(8pi)*2*pi*r*h*dr)/sum(2*pi*r*h*dr) )
            # Cancelling constant terms
            # Beavg = sqrt( sum( B**2*r*h*dr)/sum(r*h*dr) )
            quantity = (Btot**2*r*h).sum(axis=1)/(r*h).sum(axis=1)
            quantity = np.sqrt(quantity)
        else:
            #TODO
            raise NotImplementedError

    elif qname == r'B2_B2b2':
        B2  = get('Btot')*get('Btot')
        B2b2 = get('b')*get('b') + B2
        
        quantity = B2/B2b2

    elif qname == r'B2_B2b2_avg':
        if z is not None:
            B2_B2b2 = get('B2_B2b2')
            r = get('r')
            # Area weigthed average (see Bavg)
            quantity = (B2_B2b2*r).sum(axis=1)/r.sum(axis=1) #Sums for each galaxy
        else:
            #TODO
            raise NotImplementedError

    elif qname == r'Bmax_Beq':
        Bmax = get('Bmax')
        Beq = get('Beq_at_Bmax')
        quantity = Bmax/Beq

    elif qname == 'h2_at_Bmax':
        quantity = get('h_at_Bmax')*get('h_at_Bmax')

    elif qname == 'OmegaS_at_Bmax':
        quantity = get('Omega_at_Bmax')*get('Shear_at_Bmax')

    elif 'PI_I' in qname:
        match = re.match(r'PI_I(_.+cm)', qname)
        q = match.group(1)
        quantity = get('PI{}'.format(q))/get('I{}'.format(q))

    elif '_at_Bmax' in qname:
        # Extracts the name of the quantity
        match = re.match(r'(.+)_at_Bmax', qname)
        q = match.group(1)
        # Gets the quantity in the radius where B is maximum
        quantity = __get_profile_max_position(get(q), get('Btot_max_indices'),
                                              gal_id=gal_id, z=z)

    elif '_max_indices' in qname:
        # Extracts the name of the quantity
        match = re.match(r'(.+)_max_indices', qname)
        qname = match.group(1)
        # Special case: pre-computed Bmax index, stored in Bmax_idx
        # Needs to check whether the output is present in this case
        if (qname=='Btot') and ('Bmax_idx' in mag_run._data[0]):
            idx = get('Bmax_idx',0).copy()-1
            idx[np.isnan(idx)] = 0
            quantity = idx.astype(int)
        else:
            q = get(qname)
            quantity = np.argmax(q,axis=1)

    else:
        raise ValueError(qname + ' is unknown.')

    return quantity


def __get_profile_max(quantity, gal_id=None, z=None):

    if z is not None:
        quantity = quantity.max(axis=1)
    else:
        quantity = quantity.max(axis=0)
    return quantity


def __get_profile_max_position(quantity, max_indices, gal_id=None, z=None):

    if z is not None:
        indices = max_indices + np.arange(quantity.shape[0])*quantity.shape[1]
        quantity = quantity.flatten()[indices]
    else:
        #TODO
        raise NotImplementedError
    return quantity
