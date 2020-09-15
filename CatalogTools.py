from progress.bar import ChargingBar
from progress.spinner import Spinner
from astroquery.esa.hubble import ESAHubble
import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.optimize

path_sample='/Users/ziang/Documents/RAISE2020/Program/sample_hcv.csv'

path_hcv='/Users/ziang/Documents/RAISE2020/Program/1564584578352RG-result.csv'

path_monopeak='/Users/ziang/Documents/RAISE2020/Program/hcv_monopeak.csv'

def line_converter(lst):
    '''
    Given a list of raw string data and two dictionary of column-type and index-column, return a list of the correct data type.
    '''
    
    column_type={'matchid':int,'groupid':int,'subgroup':int,'ra':float,'dec':float,'pipeline_class':int,'expert_class':int,'filter':str,'num_filters':int,'var_quality_flag':str,'filter_detection_flag':int,'num_in_lc':int,'hsc_mean_mag':float,'hcv_mean_mag':float,'mad':float,'chi2':float,'lightcurve_d':float,'lightcurve_m':float,'lightcurve_cm':float,'lightcurve_e':float,'lightcurve_i':str,'lightcurve_r':str,'ci_d':float, 'ci_v':float,'d_d':float,'d_v':float}

    column_index={0: 'matchid', 1: 'groupid', 2: 'subgroup', 3: 'ra', 4: 'dec', 5: 'pipeline_class', 6: 'expert_class', 7: 'filter', 8: 'num_filters', 9: 'var_quality_flag', 10: 'filter_detection_flag', 11: 'num_in_lc', 12: 'hsc_mean_mag', 13: 'hcv_mean_mag', 14: 'mad', 15: 'chi2', 16: 'lightcurve_d', 17: 'lightcurve_m', 18: 'lightcurve_cm', 19: 'lightcurve_e', 20: 'lightcurve_i', 21: 'lightcurve_r', 22: 'ci_d', 23: 'ci_v', 24: 'd_d', 25: 'd_v'}
    
    data=[]
    for (i,val) in enumerate(lst):
        if i==18 and val=='':
            data.append(column_type[column_index[15]](lst[15]))
        else:
            try:
                data.append(column_type[column_index[i]](val))
            except ValueError:
                print(i)
                print(val)
                print(type(lst))
                break
    return(data)

def func(x,coef,c):
    return coef*(np.log10(x))+c
    
def mono_peak(lst):
    '''
    Any equality is ruled out.
    '''
    i=0
    while i<len(lst)-1 and lst[i]>lst[i+1]:
        i+=1
    
    #Get Peak.
    
    k=i
    
    while k<len(lst)-1 and lst[k]<lst[k+1]:
        k+=1
        
    return k==(len(lst)-1)

def eq_to_galactic(ra,dec):
    '''
    Given a set of equatorial coordinates, (ra,dec) in degrees, returns the J2000.0 galactic coordinates, ()
    '''
    
    ra_0=192.8595*math.pi/180
    dec_0=27.1284*math.pi/180
    l_0=122.93*math.pi/180
    ra=ra*math.pi/180
    dec=dec*math.pi/180
    
    if ra==ra_0 and dec==dec_0:
        print('This is the galactic pole. (0,0) is returned.')
        return (0,0)
    
    l_rad=l_0-math.atan(math.cos(dec)*math.sin(ra-ra_0)/(math.sin(dec)*math.cos(dec_0)-math.cos(dec)*math.sin(dec_0)*math.cos(ra-ra_0)))
    
    b_rad=math.asin(math.sin(dec)*math.sin(dec_0)+math.cos(dec)*math.cos(dec_0)*math.cos(ra-ra_0))
    
    l=l_rad*180/math.pi
    b=b_rad*180/math.pi
    
    return (l,b)
            
class Catalog:
    def __init__(self,name):
        '''
        Constructor for Catalog class. file_path is a path to a csv file.
        '''
        self.name=name
        self.catalog={}
        
    def add_target(self,target):
        '''
        Add an existing, well defined target object in this catalog.
        '''
        self.catalog[target.matchID]=target
        
    def catalog_file(self,file_path):
        try:
            print('\nReading File.')
            file=open(file_path,'r',encoding="gbk")
            hcv_lst = file.read().split("\n")[:-1]#split file into lines
            self.header=hcv_lst.pop(0).split(',')#split lines into list of str
            length=len(hcv_lst)
        
            #sort data into correct type
            bar=ChargingBar('Modifying Data Entries ',max=length)
            for i in range(length):
                hcv_lst[i]=hcv_lst[i].split(',')
                hcv_lst[i]=line_converter(hcv_lst[i])
                bar.next()
            bar.finish()
            del bar

            targetcount=0
            #sort data by target
            bar=ChargingBar('Sorting Data into Catalog ',max=length)
            i=0
            while i<length:#travese data
                current_matchID=hcv_lst[i][0]#current matchID
                target=Target(current_matchID,hcv_lst[i][1],hcv_lst[i][2],hcv_lst[i][3],hcv_lst[i][4],hcv_lst[i][5],hcv_lst[i][6]) #create new target with this matchID

                while i<length and hcv_lst[i][0]==current_matchID: #when looking at the same matchID
                    current_filter=hcv_lst[i][7] #current filter
                    lightcurve=LightCurve(hcv_lst[i][0],hcv_lst[i][7],hcv_lst[i][9],hcv_lst[i][10],hcv_lst[i][11],hcv_lst[i][12],hcv_lst[i][13],hcv_lst[i][14],hcv_lst[i][15]) #create lc with this filter
                    
                    while i<length and hcv_lst[i][7]==current_filter and hcv_lst[i][0]==current_matchID:
                        lightcurve.add_entry(Entry(hcv_lst[i])) #add an entry
                        i+=1 #cursor to next entry
                        bar.next()
                
                    target.add_data(lightcurve) #add this lightcurve to current target
                    del lightcurve
                    bar.next()
            
                #all data for this target is recorded
                self.catalog[current_matchID]=target #add this target to catalog
                targetcount+=1
                del target
                bar.next()
            bar.finish()
            del bar
            
            del hcv_lst
            
            print('Succesfully imported HCV file from '+file_path+' .')
            print('There are '+str(i)+' entries of '+str(targetcount)+' objects.\n')
            
            file.close()
        except IndentationError:
            print("Failed to generate HCV. Check file path.")
    
    def write_catalog(self,file_name):
        f=open(file_name,'w')
        
        f.truncate() 
        
        f.write('matchid,groupid,subgroup,ra,dec,pipeline_class,expert_class,filter,num_filters,var_quality_flag,filter_detection_flag,num_in_lc,hsc_mean_mag,hcv_mean_mag,mad,chi2,lightcurve_d,lightcurve_m,lightcurve_cm,lightcurve_e,lightcurve_i,lightcurve_r,ci_d,ci_v,d_d,d_v\n')
        
        for target in self.catalog.values():
            for lc in target.lightcurves:
                for entry in lc.entries:
                    line=str(target.matchID)+','+str(target.groupID)+','+str(target.subgroup)+','+str(target.ra)+','+str(target.dec)+','+str(target.pipeline_class)+','+str(target.expert_class)+','
                    
                    line+=lc.filter_type+','+str(len(target.lightcurves))+','+lc.var_quality_flag+','+str(lc.filter_detection_flag)+','+str(lc.num_in_lc)+','+str(lc.hsc_mean_m)+','+str(lc.hcv_mean_m)+','+str(lc.mad)+','+str(lc.chi2)+','
                    
                    line+=str(entry.lightcurve_d)+','+str(entry.lightcurve_m)+','+str(entry.lightcurve_cm)+','+str(entry.lightcurve_e)+','+str(entry.lightcurve_i)+','+str(entry.lightcurve_r)+','+str(entry.ci_d)+','+str(entry.ci_v)+','+str(entry.d_d)+','+str(entry.d_v)+'\n'
                    f.write(line)
                
        f.close()
        return 'Writing completed.'
                    
    def plot_time_baseline(self):
        '''
        Plot the time baseline distribution of all the targets, not lightcurves.
        See Target class for baseline definition.
        '''
        time_baseline_list=[]
        for target in self.catalog:
            time_baseline_list.append(target.get_time_baseline())
        
        if input('Plot histogram? (y/n): ')=='y':
            
            time_baseline_list_m=list(map(lambda x:x/31,time_baseline_list))
            name='Time Baseline Distribution of HCV objects'
            fig=plt.figure(num=name,figsize=(12,4),dpi=120)
            plt.hist(time_baseline_list_m,bins=np.arange(0,192,1),alpha=0.5,  histtype='stepfilled',color='steelblue',edgecolor='none',range=(0,6))
            plt.xticks(np.arange(0,192,6),fontsize=8)
            plt.yticks(np.arange(0,21000,1000),fontsize=10)
            plt.xlabel('Time baseline (Month of 31 days)',fontsize=10)
            plt.ylabel('Number of objects',fontsize=10)
            plt.title(name,fontsize=10)
            print('Time baseline histogram generated.')
            plt.show()
            if input('Save Image? (y/n): ')=='y':
                plt.savefig(name+'.png')
                print('Figure saved as "'+name+'.png"')
        return time_baseline_list
  
    def debug_1(self):
        '''
        Quick debug showing which targets don't have all the entries.
        '''
        error_lst=[]
        error_size=[]
        for ID,target in self.catalog.items():
            for lc in target.lightcurves: 
                if lc.num_in_lc !=len(lc.entries) and lc.matchID not in error_lst:
                    error_lst.append(lc.matchID)
                    error_size.append(lc.num_in_lc-len(lc.entries))
                    
        print(len(error_lst))
        print(error_lst)
        print(error_size)
        
    def entry_count(self):
        '''
        Return an integer showing the amount of entries in this catalog.
        '''
        count=0
        for ID,target in self.catalog.items():
            for lc in target.lightcurves:
                for entry in lc.entries:
                    count +=1
        return count
    
    def target_count(self):
        '''
        Return an integer showing the amount of targets in this catalog.
        '''
        count=0
        for ID,target in self.catalog.items():
            count+=1
        return count
    
    def lightcurves_count(self):
        '''
        Return an integer showing the amount of lighhtcurve in this catalog.
        '''
        count=0
        for ID,target in self.catalog.items():
            for lc in target.lightcurves:
                count+=1
        return count
    
    def filter_instr_distribution(self):
        '''
        Return a list of string recording all filter used in this catalog, and a dictionary of instrument-(list of counts) showing how many lightcurves of each filter there are taken by this instrument.
        '''
        filter_types=[]
        for key,val in self.catalog.items():
            for lc in val.lightcurves:
                n=lc.filter_type.find('_')
                #instr=lc.filter_type[:n]
                fltr=lc.filter_type[n+1:]
                    
                if fltr not in filter_types:
                    filter_types.append(fltr)
        #get a list containing the name of all filters.
        length=len(filter_types)
        
        instr_stat={}
        for key,val in self.catalog.items():
            for lc in val.lightcurves:
                n=lc.filter_type.find('_')
                instr=lc.filter_type[:n]
                fltr=lc.filter_type[n+1:]
                
                if instr not in instr_stat:
                    instr_stat[instr]=[0]*length
                    instr_stat[instr][filter_types.index(fltr)]+=1
                else:
                    instr_stat[instr][filter_types.index(fltr)]+=1
        
        return (filter_types,instr_stat)
    
    def plot_filter_distribution(self):
        '''
        Plot a bar chart, of filter distribution, with data from each instruments stacked together by same filter. Saving is optional. Does not show figure.
        '''
        (filter_types,instr_stat)=self.filter_instr_distribution()
        name='Filter and Instrument Distribution for HCV'
        
        data=[]
        label=[]
        for key,val in instr_stat.items():
            label.append(key)
            data.append(val)
            
        data=np.array(data)
        X = np.arange(data.shape[1])
        x=range(len(filter_types))
        
        fig=plt.figure(num=name,figsize=(12,4),dpi=120)
        for i in range(data.shape[0]):
            plt.bar(X,data[i],width=0.8,bottom=np.sum(data[:i],axis=0),alpha =0.8)
        plt.xticks([index for index in x],filter_types,size='small',rotation=30)
        plt.legend(label)
        plt.title(name)
        
        if input('Save image? (y/n): ')=='y':
            plt.savefig(name)
            print('Figure saved as '+name+'.png')
        
        return 'Figure generated'
         
class Target:
    def __init__(self,matchID,groupID,subgroup,ra,dec,pipeline_class,expert_class):
        '''
        Constructor for a Target class. Target class has a matchid and a list of all lightcurves as its character.
        '''
        
        self.matchID=matchID
        self.groupID=groupID
        self.subgroup=subgroup
        self.ra=ra
        self.dec=dec
        self.pipeline_class=pipeline_class
        self.expert_class=expert_class
        self.lightcurves=[]
        
    def __str__(self):
        '''
        String method for Target class. Return a string showing all key information for this entry.
        '''
        result='\n'
        result+='Target matchID: '+str(self.matchID)+'\n'
        result+='Group ID:       '+str(self.groupID)+'\n'
        result+='ra/dec:         '+str(self.ra)+', '+str(self.dec)+'\n'
        result+='Pipeline_class: '+str(self.pipeline_class)+'\n'
        result+='Expert_class:   '+str(self.expert_class)+'\n'
        result+='Num_filters:    '+str(len(self.lightcurves))+'\n'
        return result
    
    def debug_3(self):
        '''
        A quick debug showing all the image name for all its entries.
        '''
        for lc in self.lightcurves:
            for entry in lc.entries:
                print(entry.lightcurve_i)
                    
    def add_data(self,light_curve):
        '''
        Adding a light curve (Type 'LightCurve' Object) to this target.
        '''
        self.lightcurves.append(light_curve)
    
    def get_time_baseline(self):
        '''
        Return the maximum baseline for this target.
        '''
        baseline=0
        for observation in self.lightcurves:
            time_lst=[]
            for entries in observation.entries:
                time_lst.append(entries.lightcurve_d)
            
            if baseline< max(time_lst)-min(time_lst): baseline=max(time_lst)-min(time_lst)
        return baseline
    
    def get_bluest_lc(self):
        lst=[]
        for lc in self.lightcurves:
            i=lc.filter_type.find('_')
            lst.append(int(lc.filter_type[i+2:i+5]))
        
        return self.lightcurves[lst.index(min(lst))]
            
    
    def plot_lightcurves(self):
        '''
        Plot all lightcurves of this target on one figure. Save is optional.
        '''
        name='Light Curves (matchID='+str(self.matchID)+')'
        fig=plt.figure(num=name,figsize=(12,8),dpi=120)
        #xlim=(self.lightcurves[0].entries[0].lightcurve_d-50,self.lightcurves[0].entries[-1].lightcurve_d+50)
        filter_list=[]
        for lc in self.lightcurves:
            (lcx,lcy,lce)=lc.get_light_curve()
            plt.scatter(lcx, lcy,alpha=1,marker='.')
            plt.errorbar(lcx, lcy, lce,alpha=1)
            #if lc.entries[0].lightcurve_d-50 < xlim[0]: xlim[0]=lc.entries[0].lightcurve_d-50
            #if lc.entries[-1].lightcurve_d+50 > xlim[1]: xlim[1]=lc.entries[-1].lightcurve_d+50
            filter_list.append(lc.filter_type)
        plt.xlabel('MJD (days)')
        plt.ylabel('Corrected Magnitude')
        plt.title(name)
        ax=plt.gca()
        ax.xaxis.set_major_locator(plt.MultipleLocator(self.get_time_baseline()//10))
        ax.yaxis.set_major_locator(plt.MultipleLocator(0.1))
        plt.legend(filter_list)
        directory='/Users/ziang/Documents/RAISE2020/HCV/Light Curves'
        plt.savefig(directory+'/'+name+'.png')
        print('Figure saved as "'+name+'.png" to '+directory)
        return fig
                        
class LightCurve:
    def __init__(self,matchID,filter_type,var_quality_flag,filter_detection_flag,num_in_lc,hsc_mean_m,hcv_mean_m,mad,chi2):
        '''
        Constructor for LightCurve class. Instance variable contains some characteristics of this observation and the whole dataset of type 'Entry' objects.
        '''
        self.entries=[]
        self.matchID=matchID
        self.filter_type=filter_type
        self.var_quality_flag=var_quality_flag
        self.filter_detection_flag=filter_detection_flag
        self.num_in_lc=num_in_lc
        self.hsc_mean_m=hsc_mean_m
        self.hcv_mean_m=hcv_mean_m
        self.mad=mad
        self.chi2=chi2
        
    def __str__(self):
        '''
        String method for LightCurve class. Return a string showing all key information for this entry.
        '''
        result='\n'
        result+='Light Curve of matchID='+str(self.matchID)+'\n'
        result+='Filter Type:    '+self.filter_type+'\n'
        result+='Number of data: '+str(self.num_in_lc)+'\n'
        return result
    
    def get_mean_ci(self):
        i=0
        ci_sum=0
        for entry in self.entries:
            i+=1
            ci_sum+=entry.ci_v
        
        return ci_sum/i
            
    def get_baseline(self):
        time=[]
        for entry in self.entries:
            time.append(entry.lightcurve_d)
            
        return (max(time)-min(time))
               
    def add_entry(self,entry):
        '''
        A method for adding a data point to this lightcurve.
        '''
        self.entries.append(entry)
        
    def get_light_curve(self):
        '''
        Return a tuple containing the light curve data of time, magnitude and error.
        '''
        lcx=[]
        lcy=[]
        lce=[]
        for entry in self.entries:
            lcx.append(entry.lightcurve_d)
            lcy.append(entry.lightcurve_cm)
            lce.append(entry.lightcurve_e)
        
        return (lcx,lcy,lce)
    
    def well_sampling(self):
        interval=0
        for i in range(len(self.entries)-1):
            if self.entries[i+1].lightcurve_d-self.entries[i].lightcurve_d>interval:
                interval=self.entries[i+1].lightcurve_d-self.entries[i].lightcurve_d
                
        return self.get_baseline()*0.4>interval
            
    def plot_lightcurve(self):
        '''
        Plot single lightcurve for this LightCurve object. Save mendatory.
        '''
        lcx=[]
        lcy=[]
        lce=[]
        for entry in self.entries:
            lcx.append(entry.lightcurve_d)
            lcy.append(entry.lightcurve_cm)
            lce.append(entry.lightcurve_e)
        name=str(self.matchID)+' '+self.filter_type
        plt.figure(figsize=(8, 8))
        plt.scatter(lcx, lcy, color='b', alpha=0.5)
        plt.errorbar(lcx, lcy,lce, fmt='o', color='b', alpha=0.5)
        plt.xlim(min(lcx)-0.5*(max(lcx)-min(lcx)),max(lcx)+0.5*(max(lcx)-min(lcx)))
        plt.ylim(max(lcy)+0.5*(max(lcy)-min(lcy)),min(lcy)-0.5*(max(lcy)-min(lcy)))  # flip the y axis
        plt.xlabel('MJD (days)')
        plt.ylabel('Corrected Magnitude')
        plt.title('matchid='+str(self.matchID)+' filter= '+self.filter_type)
        directory='/Users/ziang/Documents/RAISE2020/HCV/hcv_TDE_model_LC'
        plt.savefig(directory+'/'+name+'.png')
        print('Figure saved as "'+name+'.png" to '+directory)
        
    def retrieve_all_images(self):
        '''
        Retrieve and save images for ALL datapoints by calling retrieve_image on all type 'Entry' objects in self.entries.
        '''
        print('Retrieve images for: matchID= '+str(self.matchID)+', filter= '+self.filter_type+'.')
        if input('There are '+str(len(self.entries))+' images to retrieve. Proceed? (y/n): ')=='y':
            for entry in self.entries:
                entry.retrieve_image()
                
    def lc_fit_exp(self):
        '''
        WIP: modify time baseline (57000+mjd seems unreal)
        WIP: modify mag (is the fall out rate in mag or log_flux or luminosity?)
        WIP: install and import packages
        '''
        x=[]
        y=[]
        for entry in self.entries:
            x.append(entry.lightcurve_d)
            y.append(entry.lightcurve_cm)
        
        popt,pcov=curve_fit(func,x,y)
        return popt[2]
    
    def lc_mono_peak(self):
        lst=[]
        for entry in self.entries:
            lst.append(entry.lightcurve_cm)
        return mono_peak(lst)
    
    def fallback_rate(self):
        y_lst=[]
        x_lst=[]
        for entry in self.entries:
            y_lst.append(entry.lightcurve_cm)
            x_lst.append(entry.lightcurve_d)
        param=scipy.optimize.curve_fit(func,x_lst,y_lst)
        return param[0]
                                          
class Entry:
    def __init__(self,line_data):
        '''
        Constructor for Entry class. Input is a list of line entry, already modified into the correct variable class.
        Instance variables includes matchID, which is essential to identifying the source this entry is for, and filter type, which sort this entry into a lightcurve, and all the data in HCV that is unique for each entry.
        '''
        self.matchID=line_data[0]
        self.filter_type=line_data[7]
        self.lightcurve_d=line_data[16]
        self.lightcurve_m=line_data[17]
        self.lightcurve_cm=line_data[18]
        self.lightcurve_e=line_data[19]
        self.lightcurve_i=line_data[20]
        self.lightcurve_r=line_data[21]
        self.ci_d=line_data[22]
        self.ci_v=line_data[23]
        self.d_d=line_data[24]
        self.d_v=line_data[25]
        
    def __str__(self):
        '''
        String method for Entry class. Return a string showing all key information for this entry.
        '''
        result+='\n'
        result+='Entry for matchID='+str(self.matchID)+'\n'
        result+='Filter and Instrument:  '+self.filter_type+'\n'
        result+='Observation date (MJD): '+str(self.lightcurve_d)+'\n'
        result+='Corrected Magnitude:    '+str(self.lightcurve_cm)+'\n'
        result+='Image File Name:        '+self.lightcurve_i+'_drz.fits\n'
        result+='Concentration Index CI: '+str(self.ci_d)+'\n'
        return result
        
    
    def retrieve_image(self):
        '''
        Retrieve and save image file shown by lightcurve_i from ESAHubble database. No return value.
        '''
        filename=self.lightcurve_i+'_drz.fits'
        print('Retrieving File '+filename)
        ESAHubble.get_artifact(filename) 
        
    def get_luminosity(self):
        L_s=3.828*(10**26)
        L=L_s*10**(0.4*(self.lightcurve_cm-4.75))
        return L
