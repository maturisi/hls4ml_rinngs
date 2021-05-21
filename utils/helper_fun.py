from include.helper_values_rings import *

def draw_PMTarray(fig, maxpmt, filemap):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    x_pmt=[]
    y_pmt=[]

    with open(filemap) as f:
        for line in f:
            line=line.strip().split()
            x_pmt.append(line[1])
            y_pmt.append(line[2])
    
    mpl.rcParams['toolbar']='None' #to disable toolbar at the bottom of figures
        #IMPORTANTE: qui va usata la mappa dove le correzioni specchi a (x,y) dei PMT siano gia state apportate
    circle_pmt=[None for _ in range(maxpmt)]
    
    for indpmt in range(len(x_pmt)):     #creating empty circles for PMTs
        if (float(x_pmt[indpmt])<2000):
            circle_pmt[indpmt]= plt.Circle((float(x_pmt[indpmt]),float(y_pmt[indpmt])), r_pmt, facecolor='none', edgecolor='black', lw=0.1)
            plt.gca().add_artist(circle_pmt[indpmt])

    return x_pmt, y_pmt


################################################

def draw_PMThits(fig, x_pmt, y_pmt, r_pmt, nhits, listch, lwidth=0.1, colorPMT='', zorderPMT=0, flag_print_chid=1):
    import matplotlib.pyplot as plt
    if colorPMT=='': #if define a color, we use that for the PMT, useful to draw only the hits used to fit a specific circle
        color_juraside='b'
        color_saleveside='k'
    else:
        color_juraside=color_saleveside=colorPMT
 

    circle_ch =[None for _ in range(maxhits_evt)]

    if (nhits!=0):
        xch=list()
        ych=list()
    else:
        print("No hits!")

    for indhits in range(nhits):        #creating filled circles for hits
        #print("chn: %d (%f,%f)" % (int(listch[indhits]), float(x_pmt[int(listch[indhits])]),float(y_pmt[int(listch[indhits])])))
        xch.append(float(x_pmt[int(listch[indhits])]))
        ych.append(float(y_pmt[int(listch[indhits])]))
        if ((xch[indhits] > 900) or (xch[indhits] > 900)): #excluding the PMT with (x,y) values outside the ARRAY, i.e. PMT not used.
            print("WARMING: PMT %d (%f, %f) outside ARRAY boundaries" % ( int(listch[indhits]), xch[indhits], ych[indhits] ))
        else:
            if int(listch[indhits])<=1024:
                circle_ch[indhits]= plt.Circle((xch[indhits],ych[indhits]), r_pmt, facecolor=color_juraside, edgecolor=color_juraside, lw=lwidth, zorder=zorderPMT)
            else:
                circle_ch[indhits]= plt.Circle((xch[indhits],ych[indhits]), r_pmt, facecolor=color_saleveside, edgecolor=color_saleveside, lw=lwidth, zorder=zorderPMT)
            plt.gca().add_patch(circle_ch[indhits])
            if (flag_print_chid): #labels with CH id near the hits circles
                plt.text(xch[indhits]+9,ych[indhits]+9, listch[indhits], fontsize=9)


################################################

def draw_Ring(fig, x_pmt, y_pmt, xc, yc, rc, listch=[], colorc='r', val_lw=1.8):
    import matplotlib.pyplot as plt
    r_pmt_reduced=4.0
    draw_PMThits(fig, x_pmt, y_pmt, r_pmt_reduced, len(listch), listch, lwidth=0.1, colorPMT='r', zorderPMT=1, flag_print_chid=0)

    circle_n= plt.Circle((xc,yc), rc, facecolor='none', edgecolor=colorc, lw=val_lw)
    plt.gca().add_patch(circle_n)
    
################################################

class Dataring:
    def __init__(self, x=0.,y=0.,r=0., tothits=0):
        self.x=x
        self.y=y
        self.r=r

################################################

def read_data_formatRECO(listfiles, all_hits_list, all_features_list, all_datarings, samples):
    import numpy as np
    hits_array=np.zeros(N_HITS_MAX)
    count_events=0
    readringdata=[]
    test_PMT=[]
    nrings=0

    for idx, file_ in enumerate(listfiles):
        with open(file_) as loadfile:
            print("opening... %s" % (file_))
            storedataring=Dataring()
            for indline,line in enumerate(loadfile):
                #print(line)
                linedata=line.split()

                if int(linedata[3])==20:
                    nrings+=1
                    
                    storedataring.x=float(linedata[5])
                    storedataring.y=float(linedata[6])
                    storedataring.r=float(linedata[7])
                    readringdata.append( storedataring )

                if int(linedata[3])==22:
                    nhits   =int(linedata[9])
                    if int(nhits)>64:
                        nhits=64
                    listch  =linedata[10:10+nhits]
                    for indh in range(len(listch)):
                        hits_array[indh]=listch[indh]

                    readringdata.append(int(nhits))
                    readringdata.append(int(nrings))
                    if nrings==0:
                        readringdata.append( storedataring ) #Se nn ci sono rings aggiungiamo struttura vuota
                    
                    all_hits_list.append(hits_array)
                    all_features_list.append(idx)
                    all_datarings.append(readringdata)

                    storedataring=Dataring()
                    nrings=0
                    hits_array=np.zeros(N_HITS_MAX)
                    count_events+=1
                if count_events==samples:
                    count_events=0
                    break

################################################

def read_data_formatNN(listfiles, all_hits_list, all_features_list, all_datarings, samples):
    import numpy as np
    debughits=0 #only to debug, it prints the hits for each event
    hits_array=np.zeros(N_HITS_MAX)
    count_events=0
    readringdata=[]
    test_PMT=[]
    
    for idx, file_ in enumerate(listfiles):
        with open(file_) as loadfile:
            print("opening... %s" % (file_))
            for indline,line in enumerate(loadfile):
                if line.split()[0].startswith('0x'):
                    #print("##### ", line.split()[0])
                    nrings=line.split()[1]
                    nhits=line.split()[3]
                    if int(nhits)>64:
                        nhits=64
                    for nh in range(int(nhits)):
                        npmt=int(line.split()[4+(3*nh)])
                        hits_array[nh]=npmt
                        if (debughits):
                            test_PMT.append(int(npmt))
                        if nh==int(nhits)-1:
                            all_hits_list.append(hits_array)
                            all_features_list.append(idx) #NOTE: data files must be ordered according to features in filelist
                            hits_array=np.zeros(N_HITS_MAX)
                            if (debughits):
                                test_PMT.sort() #sort degli hit per comodità di confronto
                                print(test_PMT)
                                test_PMT=[]
                    readringdata.append(int(nhits))
                    readringdata.append(int(nrings))
                    storedataring=Dataring()
                    #print("##--> nhits", nhits)
                    #print("##--> nrings", nrings)
                    if int(nrings) != 0:
                        for nr in range(int(nrings)):
                            ringline=next(loadfile, '').split()
                            #ringline=loadfile.next().split()[0]
                            #print("next ", ringline)
                            storedataring.x=float(ringline[2])
                            storedataring.y=float(ringline[3])
                            storedataring.r=float(ringline[4])
                            readringdata.append( storedataring )
                            storedataring=Dataring()
                    elif int(nrings)==0:
                        readringdata.append( storedataring )
                    all_datarings.append(readringdata)
                    readringdata=[] #initialize
                    count_events+=1
                    if count_events==samples:
                        count_events=0
                        break

###read_data_formatNN()

#Formato "RECO" dati usato da GL come metaformat for PATTI e dalla ricostruzione offline di NA62
#Row with "20" as third element holds ring info, 4th element number of hits used to fit, the number of PMTs is not the real ones but coming from the RECO processing before the mapping on PATTI.
#Row with "22" as third element holds events hits info corresponding to the previous number of ring(s) in the "20" line(s), 4th element is the total number of hits, the number of PMTs is relative to the RICH_remap.dat that is the PATTI mapping file.
#0 0 0 20 8 -166.951340 13.398637 163.308655 1.929874 15 1826 1868 1635 1735 1778 1872 1861 1807
#0 0 0 20 7 -213.288177 7.235085 131.247787 3.069114 15 1430 1896 1733 1889 1679 1869 1861
#0 0 0 22 15 0. 0. 0. 0. 15 515 1552 1586 523 1063 1043 1559 1546 1549 36 3 1080 1580 1564 1559

#Formato dati usato per NN (generato da acluni script come Almagest)
#0x0e074623e5 2 96.132598 14 1536 19.020000 -87.380000 1857 -304.980000 6.150000 1794 -250.980000 -25.030000 1604 -88.980000 -56.210000 1967 -430.980000 -25.030000 1770 -214.980000 -87.380000 1611 -61.980000 -196.500000 1954 -385.980000 21.740000 1743 -223.980000 -196.500000 1652 -115.980000 -196.500000 2013 -457.980000 -134.150000 1775 -250.980000 -87.380000 1692 -160.980000 -180.910000 1711 -205.980000 -102.970000
#ring: 1 -339.695645 -113.703415 127.124805 9 no_type 2013 -457.980000 -134.150000 1967 -430.980000 -25.030000 1857 -304.980000 6.150000 1794 -250.980000 -25.030000 1954 -385.980000 21.740000 1775 -250.980000 -8.380000 1743 -223.980000 -196.500000 1770 -214.980000 -87.380000 1711 -205.980000 -102.970000
#ring: 2 -95.467071 -85.697356 114.634341 4 no_type 1536 19.020000 -87.380000 1611 -61.980000 -196.500000 1692 -160.980000 -180.910000 1652 -115.980000 -196.500000
