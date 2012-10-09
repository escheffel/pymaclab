import pymaclab as pm
import pymaclab.modfiles.models as models
import pymaclab.filters as filters


def test_filters():
    rbc = pm.newMOD(models.stable.rbc1_num,mesg=True)
    rbc.modsolvers.forkleind.solve()
    # Also test here the one-element shock array at work
    rbc.modsolvers.forkleind.sim(250,('productivity'))
    
    # Use consumption and test the filtering options
    
    # First turn this into a ordinary Python list and test it
    consmat = rbc.modsolvers.forkleind.insim[1][0]
    consmat2 = consmat.reshape(consmat.shape[1],consmat.shape[0])
    consarr = consmat.__array__()
    consarr2 = consarr.flatten()
    consli = [x for x in rbc.modsolvers.forkleind.insim[1][0].__array__()[0]]
    
    datli = []
    datli.append(consmat)
    datli.append(consmat2)
    datli.append(consarr)
    datli.append(consarr2)
    datli.append(consli)

    for datar in datli:
        cons_cycle = filters.hpfilter(data=datar)[0]
        cons_cycle2 = filters.bkfilter(data=datar)
        cons_cycle3 = filters.cffilter(data=datar)[0]
