import os
import shutil

fromdir = os.path.join(
	'/','Users','mnfienen','Documents','MODELING_CENTER','NETICA_CV_GENERAL')

files = os.listdir(fromdir)
for cf in files:
	if cf.endswith('.py'):
		print u'copying over {0:s}'.format(cf)
		shutil.copy(os.path.join(fromdir,cf),cf)
		
