# mod_finder.py is used to create ω-mod files for individual species
# a FASTA file supplies the protein sequences and accession numbers to use
# Copyright © 2021 Ron Beavis

import sys
import requests
import re
import json
import os
import time
from datetime import datetime
import mysql.connector

#dict for translating common name into binomial name
SP = {	'yeast' : 'Saccharomyces cerevisiae', 
	'rat': 'Rattus norvegicus', 
	'human': 'Homo sapiens', 
	'mouse': 'Mus musculus'}

#dict for database credentials

CREDS = {	'ip':'i.i.i.i',
		'user':'uuuu',
		'pass':'pppp',
		'db': 'dddd'}


#retreive modification information for a specified accession and modification type
def getModified(_acc,_last,_type,_cursor):
	#get the GPMDB internal id number associated with the accession
	sql = "select proseqid from proseq where label=%(acc)s";
	_cursor.execute(sql,{'acc':_acc})
	vs = _cursor.fetchall()
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
	else:
		print('Bad type specification "%s"' % (_type))
	_cursor.execute(sql,{'proid':proid})
	vs = _cursor.fetchall()
	rs = dict()
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

#run through the proteins listed in _p and save information to files
def load_proteins(_p,_f,_m):
	#connect to GPMDB
	try:
		conn = mysql.connector.connect(host=CREDS['ip'],
						database=CREDS['db'],
						user=CREDS['user'],
						password=CREDS['pass') 
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
	o = open('%s_mod.xml' % (species),'w')
	js = open('%s_mod.json' % (species),'w')
	js.write('%s\n' % (json.dumps({'ENSEMBL': '104', 'species': full_species, 'GPMDB' : d, 'min obs' : _m})))
	o.write('<?xml version="1.0"?>\n')
	o.write('<bioml label="%s ENSEMBL v.104 potential modification annotation, gpmdb %s">\n' % (full_species,d))
	c = 0
	md = 0
	plen = len(_p)
	start = time.time()
	delta = time.time()
	#get modification information
	for l in _p:
		c += 1
		if c % 1000 == 0:
			dt = time.time()-delta
			tt = time.time()-start
			rate = tt/c
			rt = rate*(plen - c)
			print('%i/%i: %i (%.1f%%) delta: %.1f s, current: %.3f m, remaining: %.3f m' 
				% (c,plen,md,100.0*float(c)/float(plen),dt,tt/60,rt/60))
			delta = time.time()
		session = requests.session()
		seq = proteins[l]

		values = {'acetyl':{},'phosphoryl':{},'dimethyl':{},'GGyl':{}}
		for v in values:
			values[v] = getModified(l,len(seq),v,curses)
		min_obs = _m
		mod_types = set()
		pmods = {'accession' : l,'acetyl' : [],'phosphoryl' : [], 'dimethyl' : [], 'GGyl' : [], 'seq' : seq}
		for a in range(1,len(seq)+1):
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
		if len(mod_types):
			mods = list()
			if 'aK' in mod_types:
				mods.append('42.010565@]K')
			if 'pS' in mod_types:
				mods.append('79.966331@S')
			if 'pT' in mod_types:
				mods.append('79.966331@T')
			if 'pY' in mod_types:
				mods.append('79.966331@Y')
			if 'dR' in mod_types:
				mods.append('28.031300@R')
			line = '<protein label="%s" pmods="%s" />' % (l,','.join(mods))
			o.write('%s\n' % (line))
			md += 1
		j = json.dumps(pmods)
		js.write("%s\n" % (j))
	#disconnect from GPMDB
	curses.close()
	conn.close()
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
	f = open(sys.argv[1],'r')
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

	print('FASTA: %s\nmin obs: %i\n' % (sys.argv[1],mobs))
	load_proteins(proteins,sys.argv[1],mobs)

if __name__ == "__main__":
    main()


