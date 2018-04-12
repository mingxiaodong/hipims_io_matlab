import os, stat, re, sys, time
import struct, array
import numpy as np
import numpy.ma as ma
import matplotlib.mlab as mlab

# Mini-module written by Rob Thompson. Last update May 2014.

def readcomposite(date='20151204',timestamp='1630'):
    #inputs are date and time as strings.
    #outputs are R (rainfall rate) x and y (UK grid positions)
    
    year=date[0:4]
    month=date[4:6]
    day=date[6:8]
    
    started_at=os.getcwd()
    
    #os.chdir('/export/nuage/disk2/tmp/composite/raw/')
    #os.system('wget --no-passive-ftp ftp://ftp.ceda.ac.uk/badc/ukmo-nimrod/data/composite/uk-1km/'+year+'/metoffice-c-band-rain-radar_uk_'+date+'_1km-composite.dat.gz.tar')
    #os.system('tar xvf metoffice-c-band-rain-radar_uk_'+date+'_1km-composite.dat.gz.tar')
    #os.chdir(started_at)


    #pathed_file = '/export/nuage/disk2/tmp/composite/raw/metoffice-c-band-rain-radar_uk_'+date+timestamp+'_1km-composite.dat'
    #pathed_file = '/home/rob/Documents/Work/data/composites/metoffice-c-band-rain-radar_uk_'+date+timestamp+'_1km-composite.dat'
    print("date="+date)
    
    print("timestamp="+timestamp)
    pathed_file = 'G:\\Data\\Cumbria\\UK_Rainfall_Radar\\'+year+'\\metoffice-c-band-rain-radar_uk_'+date+timestamp+'_1km-composite.dat.gz'
    
    #print "Input file is", pathed_file
    if (os.path.isfile(pathed_file+'.gz')):
        os.system('gunzip '+pathed_file+'.gz')        

    file_id = open(pathed_file,"rb")
    record_length, = struct.unpack(">l", file_id.read(4))
    if record_length != 512: raise "Unexpected record length", record_length

    gen_ints = array.array("h")
    gen_reals = array.array("f")
    spec_reals = array.array("f")
    characters = array.array("c")
    spec_ints = array.array("h")

    gen_ints.read(file_id, 31)
    gen_ints.byteswap()
    #print "\nDate %4.4d%2.2d%2.2d Time %2.2d:%2.2d Grid %d x %d" %(gen_ints[0], gen_ints[1], gen_ints[2], gen_ints[3], gen_ints[4], gen_ints[15], gen_ints[16])

    gen_reals.read(file_id, 28)
    gen_reals.byteswap()
    #print "start northing %.1f, row interval %.1f, start easting %.1f, column interval %.1f\n"  %(gen_reals[2], gen_reals[3], gen_reals[4], gen_reals[5])

    spec_reals.read(file_id, 45)
    spec_reals.byteswap()
    characters.read(file_id, 56)
    spec_ints.read(file_id, 51)
    spec_ints.byteswap()

    record_length, = struct.unpack(">l", file_id.read(4))
    if record_length != 512: raise "Unexpected record length", record_length

    #for i in range(len(gen_ints)): print i+1, gen_ints[i]
    #for i in range(len(gen_reals)): print i+32, gen_reals[i]
    chars = characters.tostring()
    #print "Units are", chars[0:8]
    #print "Data source is", chars[8:32]
    #print "Parameter is", chars[32:55]
    #for i in range(gen_ints[22]): print i+108, spec_ints[i]

    #Read the Data
    array_size = gen_ints[15] * gen_ints[16]

    record_length, = struct.unpack(">l", file_id.read(4))
    if record_length != array_size * 2: raise "Unexpected record length", record_length

    data = array.array("h")
    try:
        data.read(file_id, array_size)
        record_length, = struct.unpack(">l", file_id.read(4))
        if record_length != array_size * 2: raise "Unexpected record length", record_length
        data.byteswap()
        #print "First 100 values are", data[:100]
        #print "Last 100 values are", data[-100:]
    except:
        print "Read failed"


    file_id.close()

    y=np.arange(gen_reals[4],gen_reals[4]+(gen_reals[5]*gen_ints[15]),gen_reals[5])/1000
    x=np.flipud(np.arange(gen_reals[2],gen_reals[2]-(gen_reals[3]*gen_ints[16]),-gen_reals[3])/1000)

    #x=[rl_datsp_hd(2)+500:rl_gen_hd(6):rl_datsp_hd(4)-500]./1000;
    #y=[rl_datsp_hd(5)+500:rl_gen_hd(4):rl_datsp_hd(1)-500]./1000;
    
    R=np.flipud(np.reshape(data,(gen_ints[15],gen_ints[16])))/32.0
    R = ma.masked_where(R<0,R)
    #Rm[Rm==0]=-10
    return R,x,y

def local_area(R,x,y,X=439.380,Y=138.564,ext=25):
    iy=mlab.find(np.abs(x-X)==min(np.abs(x-X)))
    ix=mlab.find(np.abs(y-Y)==min(np.abs(y-Y)))
    
    lR=R[ix-ext:ix+ext+1,iy-ext:iy+ext+1]
    lx=x[ix-ext:ix+ext+1]
    ly=y[iy-ext:iy+ext+1]
    
    return lR,lx,ly
    
