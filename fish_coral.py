import csv
import fractions
import os
import rasterio

from copy import deepcopy
from bs4 import BeautifulSoup
from igraph import Graph,plot
from math import nan,factorial
from matplotlib.colors import to_rgba
from mechanize import Browser
from numpy import arange,array,corrcoef,inf,linspace,ma,mean,ones,percentile,trapz,where,zeros
from random import random,randrange,sample,shuffle
from scipy.special import comb
from sklearn.ensemble import RandomForestClassifier
from time import time


def comb_(n, k):
	if k>n:
		return 0
	else:
		return factorial(n) // factorial(k) // factorial(n - k)


def cooc(j,N,N1,N2):
	a=comb_(N1,j)
	b = comb_(N-N1,N2-j)
	c = comb_(N,N2)
	return(float(fractions.Fraction(a*b,c)))



####
#RANDOM FOREST TO IDENTIFY PREY-PREDATOR PAIRS

###get eco from fishbase
fam_codes = [i[0] for i in csv.reader(open('family_codes.csv','r'))]
mech = Browser()
mech.set_handle_robots(False)

out=open('ecomat.csv','wb')
for z in fam_codes:
	url = "http://www.fishbase.org/report/KeyFactsMatrixList.php?famcode="+z
	page = mech.open(url)
	html = page.read()
	test=BeautifulSoup(html)
	span=test.findAll("span")
	if len(span)>0:
		soup = BeautifulSoup(html)
		table = soup.find("table", border=1)
		for row in table.findAll('tr')[1:]:
			col = row.findAll('td')
			sp=col[0].findAll('a')[0].renderContents().strip()
			res = [col[i].renderContents().strip().replace(b',',b'') for i in list(range(1,13))+[18]]
			out.write(b','.join([sp]+res)+b'\n')
	else:
		print ('no!!!'+z)
	print(z)

out.close()

###ecopath for depth
eco_codes = [i[0] for i in csv.reader(open('ecosystem_codes.csv','r'))]
for i in arange(0,len(eco_codes),100):
	out = open('codes_'+str(i)+'.csv','a')
	for j in eco_codes[i:i+100]:
		out.write(j[0]+'\n')
	out.close()

not_found = []
pre = [ 'https://www.fishbase.cn/',
	   'https://www.fishbase.us/',
	   'http://fishbase.mnhn.fr/',
	   'http://www.fishbase.se/',
	   'https://www.fishbase.ca/',
	   'https://www.fishbase.de/'
	  ]

out=open('drange.csv','wb')
sc = 0
for z in eco_codes:
	sc+=1
	att = 0
	while att<6:
		mech = Browser()
		mech.set_handle_robots(False)
		try:
			url = pre[att]+"Ecopath/EcoModelList.php?ve_code="+z+"&tag=1&ht=&ft="
			page = mech.open(url)
			html = page.read()
			test=BeautifulSoup(html)
			soup = BeautifulSoup(html)
			n = [i.renderContents().strip() for i in soup.findAll("div")]
			n = int([i for i in n if b"n=" in i][0].split(b'=')[1])
			if n>0:
				table = soup.find("table", border=0)
				for row in table.findAll('tr')[1:]:
					col = row.findAll('td')
					hab,fam = col[0].findAll('font')[0].renderContents().strip().replace(b'\xc2\xa0',b' ').replace(b'   ',b' ').split()
					sp=col[0].findAll('a')[0].renderContents().strip()
					res = [col[i].renderContents().strip() for i in [1,2]]
					dr = col[8].findAll('font')[0].renderContents().strip()
					out.write(b','.join([sp,hab,fam,dr]+res)+b'\n')
				print(z,len(not_found),len(eco_codes)-sc)
				done = 'yes'
				att = 10
			else:
				 att+=1
		except:
			att+=1
 	if att == 6:
		not_found.append(z)


out.close()

###############################
###download the latest interaction file from GLOBI and extract the csv file:
#https://zenodo.org/record/3950590/files/interactions.csv.gz

a = csv.reader(open('interactions.csv','r'), delimiter=',')
head = next(a)
globi = []
for row in a:
	itype, s_sp, t_sp = row[32],row[7],row[41]
	if itype in ['preysOn', 'eats'] and s_sp!='' and t_sp!='':
		#0_s_sp,1_s_gen,2_s_fam,3_s_class,4_s_phylum,5_s_kingdom,6_t_sp,7_t_gen,8_t_fam,9_t_class,10_t_phylum,11_t_kingdom,12_lat,13_lon,14_date
		row = [row[i] for i in [7,9,11,15,17,19,41,43,45,49,51,53,66,67,70]]
		globi.append(row)
	#print (len(data))


del(a)


drange = [i for i in csv.reader(open('drange.csv','r')) if i!=['finished']]
eco = [i for i in csv.reader(open('ecomat.csv','r'))]
# =============================================================================
# 0 Scientific Name
# 1 family
# 2 Max. length
# 3 L infinity
# 4 K
# 5 Years
# 6 Natural mortality
# 7 Life span
# 8 Generation time
# 9 Age at first maturity
# 10 L maturity
# 11 L max. yield
# 12 Max. weight
# 13 Trophic level
# 14 habitat
# 15 family
# 16 depth min
# 17 depth max
# 18 length
# 19 tl
# =============================================================================
dr_dict = dict([[i[0],i[1:]] for i in drange])
data_ok = []
for i in eco:
	dr = dr_dict.get(i[0],'None')
	# ['benthopelagic', 'Acestrorhynchidae', '', '30', '4.20']
	if dr!='None':
		dep = dr[2].replace(' m','').split('-')
		if dep == ['']:
			mind,maxd = '',''
		elif len(dep)==1:
			mind,maxd = float(dep[0])*0.9,float(dep[0])*1.1
		else:
			mind,maxd = dep

		data_ok.append(i+[dr[0],dr[1],mind,maxd,dr[3],dr[4]])
	else:
		data_ok.append(i+['','','','','',''])



for i in range(len(data_ok)):
	tl = data_ok[i][13]
	if tl!='':
		tl_ok = tl.replace('\xa0',' ').replace('+/-',' ').replace('s.e.',' ').split()[0]
		data_ok[i][13] = tl_ok


#fishbase lengths 'TL','WD','SL','NG','OT','DL'
all_hab = sorted(list(set([i[14].replace('\xa0','') for i in data_ok])-set([''])))
all_fam = sorted(list(set([i[1].replace('\xa0','') for i in data_ok])-set([''])))

#TL = 1.15 x SL
#TL = 1.1 x FL
data = []
for i in data_ok:
	sp = i[0].replace('\xa0','')
	fam_ = i[1].replace('\xa0','')
	fam = [1 if fam_==j else 0 for j in all_fam]
	l = i[2].replace('\xa0','')
	if 'SL' in l:
		fac = 1.15
	if 'FL' in l:
		fac = 1.1
	else:
		fac = 1
	l = float(l.replace('ot','').replace('Tl','').replace('TL','').replace('WD','').replace('SL','').replace('FL','').replace('NG','').replace('OT','').replace('DL','').replace(' ',''))*fac
	tl = i[13]
	if tl == '' or float(tl)==0:
		tl = i[19].replace('\xa0','') ###for some species, 0 means missing;check
	if tl!='' and float(tl) == 0:
		tl = ''
	hab_ = i[14].replace('\xa0','')
	hab = [1 if hab_==j else 0 for j in all_hab]
	dep_min = i[16]
	dep_max = i[17]
	len_check = i[18].replace('\xa0','')
	data.append([sp,l,tl,dep_min,dep_max]+hab)#+fam)

data_comp = [i for i in data if '' not in i]
data_comp = [[i[0]]+list(map(float,i[1:])) for i in data_comp]
eco_dict = dict([[i[0],i[1:]] for i in data_comp])
spp_ok = set([i[0] for i in data_comp])
globi_ok = [[i[0],i[6]] for i in globi if (len(set([i[0],i[6]])&spp_ok)==2)]
globi_ok = set([tuple(i) for i in globi_ok])
globi_spp = []
for i in globi_ok:
	globi_spp+=list(i)


globi_spp = set(globi_spp)
####make dataset for random forest
yx = [[1]+eco_dict[i[0]]+eco_dict[i[1]] for i in globi_ok]
pseudo_abs = []
while len(pseudo_abs)<len(yx):
	sp1,sp2 = sample(globi_spp,2)###sp1 predator,sp2 prey; keep in mind when evaluating prey pred in reef fish
	if random()<1.95:
		if tuple([sp1,sp2]) not in globi_ok and ((eco_dict[sp1][0]<eco_dict[sp2][0]*0.7) or (eco_dict[sp1][1]<eco_dict[sp2][1]*0.7)):# or (eco_dict[sp1][2]*0.5>eco_dict[sp2][3]) or (eco_dict[sp1][3]<eco_dict[sp2][2]*0.5)):
			pseudo_abs.append([0]+eco_dict[sp1]+eco_dict[sp2])
		else:
			print (sp1,sp2)
	else:
		if tuple([sp1,sp2]) not in globi_ok:
			pseudo_abs.append([0]+eco_dict[sp1]+eco_dict[sp2])
		else:
			pass#print (sp1,sp2)




yx+=pseudo_abs
shuffle(yx)
split = int(round(len(yx)*0.5))
train = yx[:split]
test = yx[split:]
train_y = [i[0] for i in train]
train_x = [i[1:] for i in train]
test_y = [i[0] for i in test]
test_x = [i[1:] for i in test]

rf = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
               max_depth=1000, max_features='auto', max_leaf_nodes=None,
               min_impurity_decrease=0.0, min_impurity_split=None,
               min_samples_leaf=1, min_samples_split=2,
               min_weight_fraction_leaf=0.0, n_estimators=1000, n_jobs=4,
               oob_score=True, random_state=0, verbose=0, warm_start=False)


rf.fit(train_x,train_y)
rf.oob_score_
pred = rf.predict_proba(test_x)
corrcoef([i[1] for i in pred],test_y)

#####after testing, calibrate the final rf on all data before
####using it to predict trophic interactions
train_y = [i[0] for i in yx]
train_x = [i[1:] for i in yx]
rf.fit(train_x,train_y)
rf.oob_score_

######################get potential fish-fish interactions from the dataset
reef_fish = sorted([i[:-4].replace('_',' ') for i in os.listdir('./ranges_fish') if i[-4:]=='.tif'])
pisc = [i[0] for i in csv.reader(open('piscivores.csv','r')) if i[0] in reef_fish]
prey = reef_fish
pisc_eco = [[i,eco_dict.get(i,None)] for i in pisc]
pisc_eco = [i for i in pisc_eco if i[1]!=None]
prey_eco = [[i,eco_dict.get(i,None)] for i in prey]
prey_eco = [i for i in prey_eco if i[1]!=None]
pot_int,pot_int_spp = [],[]
for i in pisc_eco:
	for j in prey_eco:
		if i[0]!=j[0]:
			pot_int.append(i[1]+j[1])
			pot_int_spp.append([j[0],i[0]]) ###note order: prey-->pred!!!


pred = rf.predict_proba(pot_int)
rf_prey_pred = [pot_int_spp[i]+[pred[i][1]] for i in range(len(pot_int_spp)) if pred[i][1]>0.9]


#########IDENTIFY POTENTIAL SUITABLE HABITAT FOR FISH AND CORALS
#USE MEAN TEMPERATURE - SALINITY
if not os.path.exists('./pot_ranges_fish'):
	os.makedirs('./pot_ranges_fish')

if not os.path.exists('./pot_ranges_coral'):
	os.makedirs('./pot_ranges_coral')


env = ['mean_Present.Surface.Temperature.Mean.tif','mean_Present.Surface.Salinity.Mean.tif',
'mean_bo_ph.asc','mean_Present.Surface.Primary.productivity.Mean.tif']
meow = [i for i in os.listdir('./meow_tifs') if i[-4:]=='.tif'] #rasterized regions from MEOW, see MEOW folder


env_lays = [rasterio.open('./env_rasters/'+i).read(1) for i in env]
meow_lays = [rasterio.open('./meow_tifs/'+i).read(1) for i in meow]
reef_lay = rasterio.open('reef.tif')

meta = reef_lay.meta.copy()
meta.update(compress='lzw')
reef = reef_lay.read(1)

fish_lays = sorted([i for i in os.listdir('./ranges_fish') if i[-4:]=='.tif'])
fish_spp = [i[:-4].replace('_',' ') for i in fish_lays]
for f in fish_lays:
	f1 = rasterio.open('./ranges_fish/'+f).read(1)
	if f1.sum()>0:
		occ = where(f1>0)
		pot = ones(reef.shape)*reef
		for el in env_lays:
			ov = [i for i in el[occ] if i!=nan]
			bc_min = min(ov)
			bc_max = max(ov)
			if bc_min == nan or bc_max == nan:
				print (bc_min,bc_max)
			pot[where(el<bc_min)] = 0
			pot[where(el>bc_max)] = 0
		for reg in meow_lays:
			if (f1*reg).sum() == 0:
				pot[where(reg>0)] = 0
		out=rasterio.open('./pot_ranges_fish/'+f, 'w', **meta)
		out.write(pot.astype(rasterio.uint8),1)
		out.close()


fish_op = []
for f in fish_lays:
	sp_name = f[:-4].replace('_',' ')
	f1 = rasterio.open('./ranges_fish/'+f).read(1)
	f1_pot = rasterio.open('./pot_ranges_fish/'+f).read(1)
	n = f1.sum()
	n_pot = f1_pot.sum()
	occs = where(f1>0)
	occs = ['_'.join(map(str,[occs[0][i],occs[1][i]])) for i in range(n)]
	occs_pot = where(f1_pot>0)
	occs_pot = ['_'.join(map(str,[occs_pot[0][i],occs_pot[1][i]])) for i in range(n_pot)]
	if sp_name in pisc:
		occs_null = [set(sample(occs_pot,len(occs))) for i in range(100)]
	else:
		occs_null = []
	fish_op.append([sp_name,[set(occs),occs_null]])
	print (f)



fish_op = dict(fish_op)
cooc_fish = []
sc = 0
t0 = time()
for i in rf_prey_pred:
	sp1,sp2 = fish_op[i[0]],fish_op[i[1]]
	obs_cooc = float(len(sp1[0]&sp2[0]))
	if obs_cooc>0:
		obs_null = [len(sp1[0]&sp2[1][i]) for i in range(100)]
		p = len([i for i in obs_null if i>=obs_cooc])/100.0
		if p<0.05:
			cooc_sc = (obs_cooc-mean(obs_null))/obs_cooc
			cooc_fish.append(i+[cooc_sc]) #pred prey?
			t1 = time()
			ext_t = ((t1-t0)/sc)*(len(rf_prey_pred)-sc)/3600
			print (ext_t,len(rf_prey_pred)-sc,cooc_fish[-1])
	sc+=1


###fish coral cooc
#identify coral potential ranges
coral_lays = sorted([i for i in os.listdir('./coral_ranges') if i[-4:]=='.tif'])
coral_lays_ok = []
for f in coral_lays:
	f1 = rasterio.open('./coral_ranges/'+f).read(1)
	f1 = f1*reef
	occ = where(f1>0)
	#f1[where(f1==0)] = -9999 #additional check; ok
	if f1.sum()>0:
		pot = ones(reef.shape)*reef
		for el in env_lays:
			ov = [i for i in el[occ] if i!=nan]
			bc_min = min(ov)
			bc_max = max(ov)
			pot[where(el<bc_min)] = 0
			pot[where(el>bc_max)] = 0
		for reg in meow_lays:
			if (f1*reg).sum() == 0:
				pot[where(reg>0)] = 0
		out=rasterio.open('./pot_ranges_coral/'+f, 'w', **meta)
		out.write(pot.astype(rasterio.uint8),1)
		out.close()
		coral_lays_ok.append(f)


coral_op = []
for f in coral_lays_ok:
	sp_name = f[:-4].replace('_',' ')
	f1 = rasterio.open('./coral_ranges/'+f).read(1)
	f1 = f1*reef
	f1_pot = rasterio.open('./pot_ranges_coral/'+f).read(1)
	n = f1.sum()
	n_pot = f1_pot.sum()
	occs = where(f1>0)
	occs = ['_'.join(map(str,[occs[0][i],occs[1][i]])) for i in range(n)]
	occs_pot = where(f1_pot>0)
	occs_pot = ['_'.join(map(str,[occs_pot[0][i],occs_pot[1][i]])) for i in range(n_pot)]
	occs_null = [set(sample(occs_pot,len(occs))) for i in range(100)]
	coral_op.append([sp_name,[set(occs),occs_null]])
	print (f)


coral_op = dict(coral_op)
coral_ass_f = [i for i in csv.reader(open('coral_associated.csv','r'))]
coral_ass = []
for i in coral_ass_f:
	if i[0] in fish_spp:
		w = 0
		if i[2] == 'Linked':
			w = 0.5
		if i[1]=='FAC':
			w = 0.75
		if i[1]=='OBLIG':
			w = 1
		if w!=0:
			coral_ass.append([i[0],w])



pot_fish_coral = []
for cor_gen in coral_lays_ok:
	for fish_sp,w in coral_ass:
		pot_fish_coral.append([cor_gen[:-4].replace('_',' '),fish_sp,w])




cooc_coral = []
sc = 1
t0 = time()
for i in pot_fish_coral:
	sp1,sp2 = coral_op[i[0]],fish_op[i[1]]
	obs_cooc = float(len(sp1[0]&sp2[0]))
	if obs_cooc>0:
		obs_null = [len(sp1[1][i]&sp2[0]) for i in range(100)]
		p = len([i for i in obs_null if i>=obs_cooc])/100.0
		if p<0.05:
			cooc_sc = (obs_cooc-mean(obs_null))/obs_cooc
			cooc_coral.append(i+[cooc_sc]) #pred prey?
			t1 = time()
			ext_t = ((t1-t0)/sc)*(len(pot_fish_coral)-sc)/3600
			print (ext_t,len(pot_fish_coral)-sc,cooc_coral[-1])
	sc+=1



###combine fish and coral in global_network
net = [i[:2]+[mean(i[2:])] for i in cooc_fish+cooc_coral]

net = [i for i in csv.reader(open('network.csv','r'))]
spp_in_net = []
for i in net:
	spp_in_net+=i[:2]


spp_in_net = set(spp_in_net)
###list of species per localities
reef_lay = rasterio.open('reef.tif')
meta = reef_lay.meta.copy()
meta.update(compress='lzw')
reef = reef_lay.read(1)

lat_r,lon_r = reef.shape
reef_locs = [[[] for i in range(lon_r)] for j in range(lat_r)]
for sp in spp_in_net:
	f = sp.replace(' ','_')+'.tif'
	if len(sp.split())==1:
		occ = rasterio.open('./coral_ranges/'+f).read(1)
	else:
		occ = rasterio.open('./ranges_fish/'+f).read(1)
	occ = occ*reef
	occ = where(occ>0)
	for o in range(len(occ[0])):
		reef_locs[occ[0][o]][occ[1][o]].append(sp)


div_mat_f = zeros(reef.shape)
div_mat_c = zeros(reef.shape)
for i in range(len(reef_locs)):
	for j in range(len(reef_locs[0])):
		div_mat_f[i][j] = len([k for k in reef_locs[i][j] if len(k.split())==2])
		div_mat_c[i][j] = len([k for k in reef_locs[i][j] if len(k.split())==1])



meta['dtype']='float64'
out=rasterio.open('div_mat_f.tif', 'w', **meta)
out.write(div_mat_f,1)
out.close()

out=rasterio.open('div_mat_c.tif', 'w', **meta)
out.write(div_mat_c,1)
out.close()

out=rasterio.open('div_mat.tif', 'w', **meta)
out.write(div_mat_f+div_mat_c,1)
out.close()

out = open('network.csv','w')
for i in net:
	out.write(','.join(map(str,i))+'\n')


out.close()

#############################################
####compute fish-coral dependency
net = [i for i in csv.reader(open('network.csv','r'))]
net = [i[:2]+[float(i[2])] for i in net]
net = Graph.TupleList(net,directed=True,weights=True)

reefs=where(reef>0)
res_dep = zeros(reef.shape)

out = open('results_path.csv','w')
out.close()
for reef_n in sample(range(len(reefs[0])),len(reefs[0])):
	pool = reef_locs[reefs[0][reef_n]][reefs[1][reef_n]]
	g = net.subgraph(pool)
	if len(g.es)>0:
		names = g.vs['name']
		coral_ids = [i for i in range(len(names)) if len(names[i].split())==1]
		if len(coral_ids)>0:
			fish0 = float(len([i for i in names if len(i.split())==2]))
			fish_ids = [i for i in range(len(names)) if len(names[i].split())==2]
			sss = g.shortest_paths_dijkstra(source=fish_ids, target=coral_ids, weights=None, mode='IN')
			sss_min = [min(i) for i in sss]
			net_dep = len([i for i in sss_min if i!=inf])/fish0
			cor_dep = []
			for i in g.es:
				n1,n2 = i.tuple
				if n1 in coral_ids:
					cor_dep.append(n2)
			cor_dep = len(set(cor_dep))	/fish0
			out = open('results_path.csv','a')
			out.write(','.join(map(str,[reefs[0][reef_n],reefs[1][reef_n],net_dep,cor_dep]))+'\n')
			out.close()
			print(reef_n,net_dep,cor_dep)
			res_dep[reefs[0][reef_n]][reefs[1][reef_n]] = net_dep



meta['dtype']='float64'
out=rasterio.open('net_dep.tif', 'w', **meta)
out.write(res_dep,1)
out.close()

###########################################
###################################check vulnerability to climate, habitat and fishing
#temperature tolerance
from numpy import nanmin,nanmax

env = ['mean_Present.Surface.Temperature.Mean.tif','mean_Present.Surface.Salinity.Mean.tif',
'mean_bo_ph.asc','mean_Present.Surface.Primary.productivity.Mean.tif']

env_lays = [rasterio.open('./env_rasters/'+i).read(1) for i in env]
reef_lay = rasterio.open('reef.tif')

meta = reef_lay.meta.copy()
meta.update(compress='lzw')
reef = reef_lay.read(1)

fish_lays = sorted([i for i in os.listdir('./ranges_fish') if i[-4:]=='.tif'])
fish_spp = [i[:-4].replace('_',' ') for i in fish_lays]
t_range = []
for f in fish_lays:
	f1 = rasterio.open('./ranges_fish/'+f).read(1)
	sp = f[:-4].replace('_',' ')
	if f1.sum()>0:
		occ = where(f1>0)
		pot = ones(reef.shape)*reef
		row = []
		for el in env_lays:
			ov = [i for i in el[occ] if i!=nan]
			bc_min = nanmin(ov)
			bc_max = nanmax(ov)
			row.append(bc_max-bc_min)
		t_range.append([sp,row])
	print (sp)


trange_dict = dict(t_range)
drange = [i for i in csv.reader(open('drange.csv','r')) if i!=['finished']]
eco = [i for i in csv.reader(open('ecomat.csv','r'))]
vuln = [i for i in csv.reader(open('fish_eco.csv','r'))]
#vuln[0] = ['species', 'code', 'fresh', 'brack', 'salt', 'dmin', 'dmax', 'size', 'vuln']

vuln_dict = dict([[i[0],float(i[-1])] for i in vuln[1:]])

# =============================================================================
# 0 Scientific Name
# 1 family
# 2 Max. length
# 3 L infinity
# 4 K
# 5 Years
# 6 Natural mortality
# 7 Life span
# 8 Generation time
# 9 Age at first maturity
# 10 L maturity
# 11 L max. yield
# 12 Max. weight
# 13 Trophic level
# 14 habitat
# 15 family
# 16 depth min
# 17 depth max
# 18 length
# 19 tl
# =============================================================================
dr_dict = dict([[i[0],i[1:]] for i in drange])
data_ok = []
for i in eco:
	dr = dr_dict.get(i[0],'None')
	vu = vuln_dict.get(i[0],'None')
	# ['benthopelagic', 'Acestrorhynchidae', '', '30', '4.20']
	if dr!='None':
		dep = dr[2].replace(' m','').split('-')
		if dep == ['']:
			mind,maxd = '',''
		elif len(dep)==1:
			mind,maxd = float(dep[0])*0.9,float(dep[0])*1.1
		else:
			mind,maxd = dep
		data_ok.append(i+[dr[0],dr[1],mind,maxd,dr[3],dr[4]]+[vu])
	else:
		data_ok.append(i+['','','','','','']+[vu])




#for i in range(20):
#    print (i,len([j for j in data_ok if j[i] not in ['','\xa0']]))
for i in range(len(data_ok)):
	tl = data_ok[i][13]
	if tl!='':
		tl_ok = tl.replace('\xa0',' ').replace('+/-',' ').replace('s.e.',' ').split()[0]
		data_ok[i][13] = tl_ok


#TL = 1.15 x SL
#TL = 1.1 x FL
data = []
for i in data_ok:
	sp = i[0].replace('\xa0','')
	fam = i[1].replace('\xa0','')
	l = i[2].replace('\xa0','')
	if 'SL' in l:
		fac = 1.15
	if 'FL' in l:
		fac = 1.1
	else:
		fac = 1
	l = float(l.replace('ot','').replace('Tl','').replace('TL','').replace('WD','').replace('SL','').replace('FL','').replace('NG','').replace('OT','').replace('DL','').replace(' ',''))*fac
	tl = i[13]
	if tl == '' or float(tl)==0:
		tl = i[19].replace('\xa0','') ###for some species, 0 means missing;check
	if tl!='' and float(tl) == 0:
		tl = ''
	hab = i[14].replace('\xa0','')
	dep_min = i[16]
	dep_max = i[17]
	len_check = i[18].replace('\xa0','')
	vuln = i[-1]
	data.append([sp,l,tl,dep_min,dep_max,hab,fam,vuln])



#data_comp = [i for i in data if '' not in i]
#data_comp = [[i[0]]+list(map(float,i[1:])) for i in data_comp]
eco_dict = dict([[i[0],i[1:]] for i in data])
spp_ok = set([i[0] for i in data])


env = ['mean_Present.Surface.Temperature.Mean.tif','mean_Present.Surface.Salinity.Mean.tif',
'mean_bo_ph.asc','mean_Present.Surface.Primary.productivity.Mean.tif']


reef_lay = rasterio.open('reef.tif')

meta = reef_lay.meta.copy()
meta.update(compress='lzw')
meta.update(dtype='float32')
reef = reef_lay.read(1)

vuln_mat = zeros(reef.shape)
div_mat_v = zeros(reef.shape)
div_mat_h = zeros(reef.shape)
div_mat_t = zeros(reef.shape)
div_mat = zeros(reef.shape)
hab_mat = zeros(reef.shape)
trange_mat = zeros(reef.shape)
salrange_mat = zeros(reef.shape)
phrange_mat = zeros(reef.shape)
pprange_mat = zeros(reef.shape)



fish_lays = sorted([i for i in os.listdir('./ranges_fish') if i[-4:]=='.tif'])
fish_spp = [i[:-4].replace('_',' ') for i in fish_lays]
for f in fish_lays:
	sp = f[:-4].replace('_',' ')
	sp_eco = eco_dict.get(sp,'None')
	trange = trange_dict.get(sp,'None')
	f1 = rasterio.open('./ranges_fish/'+f).read(1)
	occ = where(f1>0)
	div_mat[occ]+=1.0
	if sp_eco!='None':
		sp_vuln = sp_eco[-1]
		if sp_vuln!='None':
			if f1.sum()>0:
				vuln_mat[occ]+=sp_vuln
				div_mat_v[occ]+=1.0
		if sp_eco[4] in ['benthopelagic','demersal','reef-associated']:
				hab_mat[occ]+=1.0
		if sp_eco[4]!='':
				div_mat_h[occ]+=1.0
	if trange!='None':
		trange_mat[occ]+=trange[0]
		salrange_mat[occ]+=trange[1]
		phrange_mat[occ]+=trange[2]
		pprange_mat[occ]+=trange[3]
		div_mat_t[occ]+=1.0
	print (sp)


vuln_mat/=div_mat_v
hab_mat/=div_mat_h
trange_mat/=div_mat_t
salrange_mat/=div_mat_t
phrange_mat/=div_mat_t
pprange_mat/=div_mat_t

vuln_mat[where(vuln_mat==nan)] = 0
hab_mat[where(hab_mat==nan)] = 0
trange_mat[where(trange_mat==nan)] = 0
salrange_mat[where(salrange_mat==nan)] = 0
phrange_mat[where(phrange_mat==nan)] = 0
pprange_mat[where(pprange_mat==nan)] = 0


out=rasterio.open('vuln_mat.tif', 'w', **meta)
out.write(vuln_mat.astype(rasterio.float32),1)
out.close()


out=rasterio.open('hab_mat.tif', 'w', **meta)
out.write(hab_mat.astype(rasterio.float32),1)
out.close()

out=rasterio.open('trange_mat.tif', 'w', **meta)
out.write(trange_mat.astype(rasterio.float32),1)
out.close()

out=rasterio.open('salrange_mat.tif', 'w', **meta)
out.write(salrange_mat.astype(rasterio.float32),1)
out.close()

out=rasterio.open('phrange_mat.tif', 'w', **meta)
out.write(phrange_mat.astype(rasterio.float32),1)
out.close()

out=rasterio.open('pprange_mat.tif', 'w', **meta)
out.write(pprange_mat.astype(rasterio.float32),1)
out.close()

out=rasterio.open('fish_div.tif', 'w', **meta)
out.write(div_mat.astype(rasterio.float32),1)
out.close()


#######explore dependency vs other factors

reef_lay = rasterio.open('reef.tif')
meta = reef_lay.meta.copy()
meta.update(compress='lzw')
reef = reef_lay.read(1)

res_dep = rasterio.open('net_dep.tif').read(1)
iso = rasterio.open('./env_rasters/isolation.tif').read(1)
reef_acc =  rasterio.open('./env_rasters/reef_acc.tif').read(1)

fish_vuln = rasterio.open('vuln_mat.tif').read(1)
fish_hab = rasterio.open('hab_mat.tif').read(1)
fish_trange = rasterio.open('trange_mat.tif').read(1)
fish_salrange = rasterio.open('salrange_mat.tif').read(1)
fish_phrange = rasterio.open('phrange_mat.tif').read(1)
fish_pprange = rasterio.open('pprange_mat.tif').read(1)



fish_div = rasterio.open('fish_div.tif').read(1)

fff = os.listdir('./ocean_impact_layers/')

glob_driv = ['oa_2013_impact.tif',
 'slr_2013_impact.tif',
 'sst_2013_impact.tif']

hum_driv = ['light_2013_impact.tif',
 'shipping_2013_impact.tif',
 'commercial_fishing_pelagic_high_bycatch_2013_impact.tif',
 'commercial_fishing_pelagic_low_bycatch_2013_impact.tif',
 'commercial_fishing_demersal_nondestructive_low_bycatch_2013_impact.tif',
 'direct_human_2013_impact.tif',
 'commercial_fishing_demersal_nondestructive_high_bycatch_2013_impact.tif',
 'artisanal_fishing_2013_impact.tif',
 'nutrient_pollution_2013_impact.tif',
 'commercial_fishing_demersal_destructive_2013_impact.tif',
 'organic_chemical_pollution_2013_impact.tif']

#h_imp =  rasterio.open('./env_rasters/others/human_impact_reef.tif').read(1)
h_imp = 'none'
for lay in hum_driv:
	if h_imp == 'none':
		h_imp = rasterio.open('./ocean_impact_layers/'+lay).read(1)
	else:
		h_imp += rasterio.open('./ocean_impact_layers/'+lay).read(1)

h_imp*=reef

glob_imp = 'none'
for lay in glob_driv:
	if glob_imp == 'none':
		glob_imp = rasterio.open('./ocean_impact_layers/'+lay).read(1)
	else:
		glob_imp += rasterio.open('./ocean_impact_layers/'+lay).read(1)


glob_imp*=reef

tot_imp = h_imp+glob_imp
ta_ssp2 =  rasterio.open('./env_rasters/ssp2_ta.tif').read(1)
ta_ssp3 =  rasterio.open('./env_rasters/ssp3_ta.tif').read(1)
ta_ssp5 =  rasterio.open('./env_rasters/ssp5_ta.tif').read(1)

bleach =  rasterio.open('./env_rasters/bleaching.tif').read(1)
oblig =  rasterio.open('./env_rasters/oblig_dependency.tif').read(1)
fac =  rasterio.open('./env_rasters/fac_dependency.tif').read(1)
linked =  rasterio.open('./env_rasters/linked_dependency.tif').read(1)

out = open('net_dep_vs_acc.csv','w')
out.write('net_dep,acc,iso,h_imp,glob_imp,bleach,ta_ssp2,ta_ssp3,ta_ssp5,oblig,fac,linked,fish_vuln,hab,trange,salrange,phrange,pprange,div\n')
for i in range(len(res_dep)):
	for j in range(len(res_dep[0])):
		if res_dep[i][j]!=0:
			out.write(','.join(map(str,[res_dep[i][j],
									   reef_acc[i][j],
									   iso[i][j],
									   h_imp[i][j],
									   glob_imp[i][j],
									   bleach[i][j],
									   ta_ssp2[i][j],
									   ta_ssp3[i][j],
									   ta_ssp5[i][j],
									   oblig[i][j],
									   fac[i][j],
									   linked[i][j],
									   fish_vuln[i][j],
									   fish_hab[i][j],
									   fish_trange[i][j],
									   fish_salrange[i][j],
									   fish_phrange[i][j],
									   fish_pprange[i][j],
									   fish_div[i][j]
									   ]))+'\n')


out.close()


reef_lay = rasterio.open('reef.tif')
meta = reef_lay.meta.copy()
meta.update(compress='lzw')
reef = reef_lay.read(1)

imp_t = percentile(h_imp[where(reef>0)],70)
dep_t = percentile(res_dep[where(reef>0)],70)
bleach_t = percentile(bleach[where(reef>0)],70)
glob_t = percentile(glob_imp[where(reef>0)],70)
tot_t = percentile(tot_imp[where(reef>0)],70)

h_imp_norm = 1.0*(h_imp>=imp_t)#(h_imp-h_imp.min())/(h_imp.max()-h_imp.min())
tot_imp_norm = 1.0*(tot_imp>=tot_t)#(h_imp-h_imp.min())/(h_imp.max()-h_imp.min())

dep_norm =1.0*(bleach>=bleach_t)*(res_dep>=dep_t) #(bleach-bleach.min())/(bleach.max()-bleach.min())
glob_imp_norm = 1.0*(glob_imp>=glob_t)#(h_imp-h_imp.min())/(h_imp.max()-h_imp.min())


# =============================================================================
tot_vuln_1 = zeros(reef.shape)
for i in range(len(reef)):
	for j in range(len(reef[0])):
 		tot_vuln_1[i][j] = max(dep_norm[i][j],tot_imp_norm[i][j])
		#tot_vuln_1[i][j] = max(glob_imp_norm[i][j],h_imp_norm[i][j])

# =============================================================================
tot_vuln = 2*tot_imp_norm+dep_norm#glob_imp_norm

h_imp_norm.sum()
glob_imp_norm.sum()
dep_norm.sum()
tot_imp_norm.sum()
tot_vuln_1.sum()

reef.sum()

meta['dtype']='float64'
out=rasterio.open('tot_vuln.tif', 'w', **meta)
out.write(tot_vuln,1)
out.close()

meta['dtype']='float32'

out=rasterio.open('human_impact_reef.tif', 'w', **meta) ##local hazards
out.write(h_imp,1)
out.close()

out=rasterio.open('glob_imp.tif', 'w', **meta)#global hazards
out.write(glob_imp,1)
out.close()


out=rasterio.open('tot_imp.tif', 'w', **meta)#total hazards
out.write(glob_imp+h_imp,1)
out.close()


#out=rasterio.open('bleaching_norm.tif', 'w', **meta)
#out.write(bleach_norm,1)
#out.close()

out=rasterio.open('net_dep_norm.tif', 'w', **meta)
out.write(dep_norm,1)
out.close()



out = open('tot_vuln.csv','w')
for i in range(len(reef)):
	for j in range(len(reef)):
		if reef[i][j]>0:
			out.write(','.join(map(str,[h_imp_norm[i][j],tot_vuln[i][j]]))+'\n')


out.close()


####explore sensitivity to percentile choice

reef_n = reef.sum()

out = open('ov_vs_plev.csv','w')
out.write('p_lev,imp_n,dep_n,ov,tot_vuln,imp_n%,dep_n%,increase\n')
for p_lev in arange(0,101,1):
	imp_t = percentile(tot_imp[where(reef>0)],p_lev)
#	glob_t = percentile(glob_imp[where(reef>0)],p_lev)
	dep_t = percentile(res_dep[where(reef>0)],p_lev)
	bleach_t = percentile(bleach[where(reef>0)],p_lev)
	h_imp_norm = 1.0*(tot_imp>=imp_t)#(h_imp-h_imp.min())/(h_imp.max()-h_imp.min())
#	glob_imp_norm = 1.0*(glob_imp>=glob_t)#(h_imp-h_imp.min())/(h_imp.max()-h_imp.min())
	dep_norm =1.0*(bleach>=bleach_t)*(res_dep>=dep_t) #(bleach-bleach.min())/(bleach.max()-bleach.min())
	tot_vuln_1 = zeros(reef.shape)
	for i in range(len(reef)):
		for j in range(len(reef[0])):
	 		tot_vuln_1[i][j] = max(dep_norm[i][j],h_imp_norm[i][j])#max(glob_imp_norm[i][j],h_imp_norm[i][j])#
	imp_n = h_imp_norm.sum()
	dep_n = dep_norm.sum()#glob_imp_norm.sum()
	ov = ((imp_n+dep_n)-tot_vuln_1.sum())/imp_n
	out.write(','.join(map(str,[round(i,2) for i in [p_lev,imp_n,dep_n,ov,tot_vuln_1.sum()/reef_n,imp_n/reef_n,dep_n/reef_n,(tot_vuln_1.sum()-imp_n)/imp_n]]))+'\n')


out.close()






