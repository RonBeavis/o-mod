# mod_finder.py is used to create ω-mod files for individual species
# a FASTA file supplies the protein sequences and accession numbers to use
# Copyright © 2021 Ron Beavis
# Note: this file is for production only and it should not be shared with
# correct databases CREDS.

import sys
import re
import json
import os
import time
import re
from datetime import datetime
import mysql.connector

#dict for translating common name into binomial name
SP = {	'yeast' : 'Saccharomyces cerevisiae', 
	'rat': 'Rattus norvegicus', 
	'human': 'Homo sapiens', 
	'mouse': 'Mus musculus'}


#in this case, anonymized and non-functional
#CREDS = {	'ip':'i.i.i.i',
#		'user':'uuuu',
#		'pass':'pppp',
#		'db': 'dddd'}

GENES = {	'semi': set(['ins','ins1','ins2']),
		'gla': set(['f2','pros1','f10','f9','f7','proc','proz','bglap']),
		'gG' : set(['uba52','rps27a','ubb','ubc']),
		'cP' : set([	'fcn1','fcn2','fcn3','fga','c1qa','c1qb','c1qc','adipoq','sftpb',
				'sftpa2','sftpa1','sftpc','colq','marco','gldn','emid1','emid2',
				'emilin1','emilin2']),
		'iY' : set(['tg'])
	}

def getDescription(_acc,_cursor):
	sql = "select description from enspmapdb.map where label=%(label)s"
	_cursor.execute(sql,{'label':_acc})
	vs = _cursor.fetchall()
	if vs:
		return vs[0][0]
	return ""

#retreive modification information for a specified accession and modification type
def getModified(_acc,_last,_type,_cursor):
	#get the GPMDB internal id number associated with the accession
	sql = "select proseqid from proseq where label=%(acc)s";
	_cursor.execute(sql,{'acc':_acc})
	vs = _cursor.fetchall()
	#if the accession number does not exist, return 0 for all residues
	if not vs:
		res = dict()
		for i in range(1,_last+1):
			res[i] = 0
		return res

	proid = vs[0][0]
	#setup the SQL statement for getting modification information
	if _type == 'phosphoryl':
		sql = "select at,freq from phos_freq where proseqid=%(proid)s";
	elif _type == 'acetyl':
		sql = "select at,freq from acetyl_freq where proseqid=%(proid)s";
	elif _type == 'dimethyl':
		sql = "select at,freq from dim_freq where proseqid=%(proid)s";
	elif _type == 'GGyl':
		sql = "select at,freq from ubi_freq where proseqid=%(proid)s";
	elif _type == 'hydroxyl':
		sql = "select at,freq from oxy_freq where proseqid=%(proid)s";
	else:
		print('Error:bad type specification "%s"' % (_type))
		exit()
	_cursor.execute(sql,{'proid':proid})
	vs = _cursor.fetchall()
	rs = dict()
	#create a dict for non-zero values
	for v in vs:
		rs[v[0]] = v[1]
	res = dict()
	#create a dict with values across the entire sequence
	for i in range(1,_last+1):
		if i in rs:
			res[i] = rs[i]
		else:
			res[i] = 0
	return res

def check_gene(_d):
	t = re.sub(r'\:p.+',r'',_d.lower())
	vs = []
	for x in GENES:
		if t in GENES[x]:
			print('\t%s, %s' % (x,_d))
			vs.append(x)
	ls = re.findall('^krt\d+$',t)
	if len(ls):
		x = 'cR'
		print('\t%s, %s' % (x,_d))
		vs.append(x)
	ls = re.findall('^col\d+a\d$',t)
	if len(ls):
		x = 'cP'
		print('\t%s, %s' % (x,_d))
		vs.append(x)
	return vs

#run through the proteins listed in _p and save information to files
def load_proteins(_p,_f,_m):
	#connect to GPMDB
	try:
		conn = mysql.connector.connect(host=CREDS['ip'],
						database=CREDS['db'],
						user=CREDS['user'],
						password=CREDS['pass']) 
	except:
		print('{Error: Could not connect to database')
		exit()
	#create a cursor
	curses = conn.cursor()
	#get species common name from fasta file
	species = re.sub(r'(.+?)_.+',r'\1',_f)
	full_species = species
	#get binomial name from SP
	if species in SP:
		full_species = SP[species]
	#compose file header information
	d = datetime.today().strftime('%Y-%m-%d')
	o = open('mods/%s_mod.xml' % (species),'w')
	js = open('json/%s_mod.json' % (species),'w')
	js.write('%s\n' % (json.dumps({'ENSEMBL': '104', 'species': full_species, 'GPMDB' : d, 'min obs' : _m})))
	o.write('<?xml version="1.0"?>\n')
	o.write('<bioml label="%s ENSEMBL v.104 potential modification annotation, gpmdb %s">\n' % (full_species,d))
	c = 0
	md = 0
	plen = len(_p)
	#initialize timers to monitor progress
	start = time.time()
	delta = time.time()
	#get modification information
	for l in _p:
		c += 1
		# print a keep-alive diagnostic message at regular intervals
		if c % 1000 == 0:
			dt = time.time()-delta
			tt = time.time()-start
			rate = tt/c
			rt = rate*(plen - c)
			print('%i/%i: %i (%.1f%%) delta: %.1f s, current: %.3f m, remaining: %.3f m' 
				% (c,plen,md,100.0*float(c)/float(plen),dt,tt/60,rt/60))
			delta = time.time()
		seq = _p[l]
		#dict for storing results, by PTM
		values = {'acetyl':{},'phosphoryl':{},'dimethyl':{},'GGyl':{},'hydroxyl':{}}
		#load values
		for v in values:
			values[v] = getModified(l,len(seq),v,curses)
		min_obs = _m
		mod_types = set()
		#initialize the dict that will be used to generate single line JSON object
		desc = getDescription(l,curses)
		pmods = {'accession' : l, 'description': desc, 'acetyl' : [],'phosphoryl' : [], 'dimethyl' : [], 'GGyl' : [], 'seq' : seq}
		#loop through all residues, updating records for each
		hydroxyl = 0
		for a in range(1,len(seq)+1):
			if a in values['hydroxyl']:
				if seq[a-1] == 'P' or seq[a-1] == 'K':
					p = values['hydroxyl'][a]
					if p >= 20:
						mod_types.add('oP')
						hydroxyl += 1
			if a in values['acetyl']:
				if seq[a-1] == 'K':
					p = values['acetyl'][a]
					if p >= min_obs:
						mod_types.add('aK')
					if p > 0:
						pmods['acetyl'].append({'r':seq[a-1],'c':a,'n':p})
			if a in values['GGyl']:
				if seq[a-1] == 'K':
					p = values['GGyl'][a]
					if p >= min_obs:
						mod_types.add('gK')
					if p > 0:
						pmods['GGyl'].append({'r':seq[a-1],'c':a,'n':p})
			if a in values['dimethyl']:
				if seq[a-1] == 'R':
					p = values['dimethyl'][a]
					if p >= min_obs:
						mod_types.add('dR')
					if p > 0:
						pmods['dimethyl'].append({'r':seq[a-1],'c':a,'n':p})
			if a in values['phosphoryl']:
				p = values['phosphoryl'][a]
				if seq[a-1] == 'S':
					if p >= min_obs:
						mod_types.add('pS')
					if p > 0:
						pmods['phosphoryl'].append({'r':seq[a-1],'c':a,'n':p})
				if seq[a-1] == 'T':
					if p >= min_obs:
						mod_types.add('pT')
					if p > 0:
						pmods['phosphoryl'].append({'r':seq[a-1],'c':a,'n':p})
				if seq[a-1] == 'Y':
					if p >= min_obs:
						mod_types.add('pY')
					if p > 0:
						pmods['phosphoryl'].append({'r':seq[a-1],'c':a,'n':p})
		vs = check_gene(desc)
		for v in vs:
			mod_types.add(v)
		if vs:
			print('\tmods, ',mod_types)
		#check to see if the XML file will be updated
		if len(mod_types):
			mods = set()
			if 'aK' in mod_types:
				mods.add('42.010565@]K')
			if 'pS' in mod_types:
				mods.add('79.966331@S')
			if 'pT' in mod_types:
				mods.add('79.966331@T')
			if 'pY' in mod_types:
				mods.add('79.966331@Y')
			if 'dR' in mod_types:
				mods.add('28.031300@]R')
			if 'cR' in mod_types:
				mods.add('0.984016@]R')
			if 'oP' in mod_types and hydroxyl > 4:
				mods.add('15.994915@K')
				mods.add('15.994915@P')
				mods.add('-1@B')
			elif 'cP' in mod_types:
				mods.add('15.994915@K')
				mods.add('15.994915@P')
				mods.add('-1@B')
			if 'gla' in mod_types:
				mods.add('43.989829@E')
			if 'iY' in mod_types:
				mods.add('125.896648@Y')
			if 'gG' in mod_types:
				mods.add('114.042927@]K')
				mods.add('-1@B')
			if 'semi' in mod_types:
				mods.add('-1@X')
			pl = ','.join(mods)
			if len(pl) > 0:
				line = '<protein label="%s" pmods="%s" />' % (l,pl)
				o.write('%s\n' % (line))
			md += 1
		#record the information in pmods as a line of JSON
		j = json.dumps(pmods)
		js.write("%s\n" % (j))
	#disconnect from GPMDB
	curses.close()
	conn.close()
	#add additional hand-curated information
	if os.path.exists('%s_extra.xml' % (species)):
		print('Adding %s_extra.xml' % (species))
		ls = [l.strip() for l in open('%s_extra.xml' % (species),'r')]
		for line in ls:
			if line.find('<') == -1:
				continue
			o.write('%s\n' % (line))
	else:
		print('No %s_extra.xml exists' % (species))

	#compose and record file finish information

	o.write('</bioml>\n')
	o.close()
	
	js.write("%s\n" % (json.dumps({'lines': c+2, 'accessions': c, 'modified': md})))
	js.close()

def main():
	#open the specified FASTA file
	base = sys.argv[1]
	f = open('fasta/%s' % (base),'r')
	proteins = dict()
	label = ''
	seq = ''
	#store protein information into proteins
	for l in f:
		l = l.strip()
		if l.find('>ENS') == 0:
			if len(label):
				proteins[label] = seq
			label = re.sub(r'\>(ENS[A-Z]*\d+).+',r'\1',l)
			seq = ''
		elif l.find('>') == 0:
			if len(label):
				proteins[label] = seq
			label = re.sub(r'\>([^ ]+) .+',r'\1',l)
			seq = ''
		else:
			seq += l
	if len(seq):
		proteins[label] = seq
	f.close()
	#retreive minimum # of observations from command line
	mobs = 5
	try:
		mobs = int(sys.argv[2])
	except:
		pass

	#start getting protein modification information

	print('FASTA: %s\nmin obs: %i\n' % (base,mobs))
	load_proteins(proteins,base,mobs)

if __name__ == "__main__":
    main()


