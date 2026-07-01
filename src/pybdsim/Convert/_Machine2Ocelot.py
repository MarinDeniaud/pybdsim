import numpy as _np
import pybdsim as _bd
import ocelot as _ocl


def Machine2Ocelot(bdsmachine, s0=0):
    """ Convert BDSIM machine to an Ocelot lattice.

       +-----------------+---------------------------------------------------------+
       | **Parameters**  | **Description**                                         |
       +-----------------+---------------------------------------------------------+
       | bdsmachine      | BDSIM machine element                                   |
       +-----------------+---------------------------------------------------------+
       | s0              | Initial s position of the line. Optional.               |
       +-----------------+---------------------------------------------------------+
    """

    def loadParam(machine, key):
        try:
            string = machine.beam[key]
        except KeyError:
            return 0.0
        if '*' in string:
            value_str, unit = string.split("*")
            if key == 'energy':
                if unit == 'eV':
                    return float(value_str)*1e-9
                elif unit == 'keV':
                    return float(value_str)*1e-6
                elif unit == 'MeV':
                    return float(value_str)*1e-3
                elif unit == 'GeV':
                    return float(value_str)
            else:
                return float(value_str)
        else:
            if key == 'energy':
                return float(string)*1e-9
            return float(string)

    # Initial Twiss parameters
    tws0 = _ocl.Twiss(beta_x=loadParam(bdsmachine, 'betx'), beta_y=loadParam(bdsmachine, 'bety'),
                      alpha_x=loadParam(bdsmachine, 'alfx'), alpha_y=loadParam(bdsmachine, 'alfy'),
                      Dx=loadParam(bdsmachine, 'dispx'), Dy=loadParam(bdsmachine, 'dispy'),
                      Dxp=loadParam(bdsmachine, 'dispxp'), Dyp=loadParam(bdsmachine, 'dispyp'),
                      mux=loadParam(bdsmachine, 'mux'), muy=loadParam(bdsmachine, 'muy'),
                      emit_x=loadParam(bdsmachine, 'emitx'), emit_y=loadParam(bdsmachine, 'emity'),
                      E=loadParam(bdsmachine, 'energy'), pp=loadParam(bdsmachine, 'sigmaE'),
                      s=s0)

    cell = []
    # unique_name_list = []

    for elem_name in bdsmachine.sequence:
        # if elem_name not in unique_name_list:
        #     unique_name_list.append(elem_name)

        element = bdsmachine.elements[elem_name]

        factor = 1
        if bdsmachine.charge == -1:
            factor = -1

        match element.category:
            case 'marker':
                cell.append(_ocl.Marker(eid=elem_name))
            case 'drift':
                cell.append(_ocl.Drift(eid=elem_name, l=element.length))
            case 'rbend':
                e1 = element.get('e1') if element.get('e1') is not None else 0
                e2 = element.get('e2') if element.get('e2') is not None else 0
                k1 = element.get('k1') if element.get('k1') is not None else 0
                cell.append(_ocl.RBend(eid=elem_name, l=element.length, angle=element.get('angle'),
                                       e1=e1, e2=e2, k1=k1))
            case 'sbend':
                e1 = element.get('e1') if element.get('e1') is not None else 0
                e2 = element.get('e2') if element.get('e2') is not None else 0
                k1 = element.get('k1') if element.get('k1') is not None else 0
                cell.append(_ocl.SBend(eid=elem_name, l=element.length, angle=element.get('angle'),
                                       e1=e1, e2=e2, k1=k1))
            case 'quadrupole':
                cell.append(_ocl.Quadrupole(eid=elem_name, l=element.length, k1=factor*element.get('k1')))
            case 'sextupole':
                cell.append(_ocl.Sextupole(eid=elem_name, l=element.length, k2=factor*element.get('k2')))
            case 'octupole':
                cell.append(_ocl.Octupole(eid=elem_name, l=element.length, k3=factor*element.get('k3')))
            case 'solenoid':
                cell.append(_ocl.Solenoid(eid=elem_name, l=element.length, k=element.get('ks')))
            case 'rfcavity':
                volt = element.length * element.get('gradient')
                cell.append(_ocl.Cavity(eid=elem_name, l=element.length, freq=element.get('freq'),
                                        v=volt, phi=element.get('phase')))
            case _:
                print('Unknown element type:', element.category)
                if element.length > 0:
                    cell.append(_ocl.Drift(eid=elem_name, l=element.length))
                else:
                    cell.append(_ocl.Marker(eid=elem_name))
        if element.get('tilt') is not None and element.get('tilt') > 0:
            cell[-1].tilt = element.get('tilt')

    return _ocl.MagneticLattice(cell, method={'global': _ocl.SecondTM}), tws0
