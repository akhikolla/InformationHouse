from __future__ import print_function
from __future__ import unicode_literals
from builtins import str
import os
import sys
import subprocess
import string
import argparse
from time import gmtime, strftime

version = 'v2.0'
print('createImageSeq {}'.format(version))
print('python {}'.format(sys.version))

python3 = sys.version_info > (3, 0)

nodate_prefix = 'nodate_'

stand_alone_testing = False
	
def conv_command(fileName, fps, hor, vert, targetDirName, start_given, stop_given):
	inFile = os.path.join(movieDir, fileName)
	outFiles = os.path.join(sequenceParentDir, targetDirName, 'image-%05d.png')
	return [libavpath,
			'-loglevel', 'quiet',
			'-i', inFile,
			'-r', str(fps),
			'-s', str(hor) + 'x' + str(vert),
			'-ss', str(start_given),
			'-t', str(stop_given - start_given),
			'-f', 'image2',
			outFiles]
def which(cmd, mode=os.F_OK | os.X_OK, path=None):
    """Given a command, mode, and a PATH string, return the path which
    conforms to the given mode on the PATH, or None if there is no such
    file.
    `mode` defaults to os.F_OK | os.X_OK. `path` defaults to the result
    of os.environ.get("PATH"), or can be overridden with a custom search
    path.
	Copied from 3.6 source referred to at https://docs.python.org/3.6/library/shutil.html,
	also works for 2.7.
    """
    # Check that a given file can be accessed with the correct mode.
    # Additionally check that `file` is not a directory, as on Windows
    # directories pass the os.access check.
    def _access_check(fn, mode):
        return (os.path.exists(fn) and os.access(fn, mode)
                and not os.path.isdir(fn))
    # If we're given a path with a directory part, look it up directly rather
    # than referring to PATH directories. This includes checking relative to the
    # current directory, e.g. ./script
    if os.path.dirname(cmd):
        if _access_check(cmd, mode):
            return cmd
        return None
    if path is None:
        path = os.environ.get('PATH', os.defpath)
    if not path:
        return None
    path = path.split(os.pathsep)
    if sys.platform == u"win32":
        # The current directory takes precedence on Windows.
        if not os.curdir in path:
            path.insert(0, os.curdir)
        # PATHEXT is necessary to check on Windows.
        pathext = os.environ.get("PATHEXT", "").split(os.pathsep)
        # See if the given file matches any of the expected path extensions.
        # This will allow us to short circuit when given "python.exe".
        # If it does match, only test that one, otherwise we have to try
        # others.
        if any(cmd.lower().endswith(ext.lower()) for ext in pathext):
            files = [cmd]
        else:
            files = [cmd + ext for ext in pathext]
    else:
        # On other platforms you don't have things like PATHEXT to tell you
        # what file suffixes are executable, so just pass on cmd as-is.
        files = [cmd]
    seen = set()
    for dir in path:
        normdir = os.path.normcase(dir)
        if not normdir in seen:
            seen.add(normdir)
            for thefile in files:
                name = os.path.join(dir, thefile)
                if _access_check(name, mode):
                    return name
    return None
parser = argparse.ArgumentParser()
parser.add_argument('-moviepath')
parser.add_argument('-imagepath')
parser.add_argument('-libavpath')
parser.add_argument('-exiftoolpath')
parser.add_argument('-x')
parser.add_argument('-y')
parser.add_argument('-fps')
parser.add_argument('-nsec')
parser.add_argument('-start')
parser.add_argument('-stop')
parser.add_argument('-ext', nargs='*')
parser.add_argument('-verbose')
parser.add_argument('-logfile')
if stand_alone_testing:
	from argparse import Namespace
	args = Namespace()
	args.moviepath = 'Movies'
	args.imagepath = 'ImageSequences'
	args.libavpath = 'avconv'
	args.exiftoolpath = 'exiftool'
	args.x = '1920'
	args.y = '1080'
	args.fps = '15'
	args.nsec = '2'
	args.start = 'NULL'
	args.stop = 'NULL'
	args.ext = ['mts', 'avi']
	args.verbose = 'TRUE'
	args.logfile = 'TRUE'
else:
	args = parser.parse_args()
if args.logfile not in ['FALSE', 'TRUE']:
	print('Unknown value %s for logfile. Converting to FALSE' % args.verbose)
	logtofile = False
else:
	logtofile = args.logfile == 'TRUE'
if logtofile:
	timeofday = gmtime()
	logfilename = os.path.join(os.getcwd(), 'trackdem_' + strftime("%Y-%m-%d_%H-%M-%S", timeofday) + '.log')
	logfile = open(logfilename, 'w')
	logfile.write('CreateImageSeq ' + version + ' log created on ' + strftime("%Y-%m-%d %H:%M:%S", timeofday) + '\n\n')
	print('Writing to log file ' + logfilename)
if logtofile:
	logfile.write('args = %s\n\n' % str(args))
# process arguments
libavpath = args.libavpath
exiftoolpath = args.exiftoolpath
if logtofile:
	logfile.write('Calling "which" with command: "%s"\n\n' % 'libavpath')
avconv = which(libavpath)
if logtofile:
	logfile.write('Calling "which" with command: "%s"\n\n' % 'exiftool')
exiftool = which(exiftoolpath)
if avconv:
	if logtofile:
		logfile.write('Using avconv found at %s\n\n' % avconv)
else:
	if logtofile:
		logfile.write('Cannot find avconv executable, exiting\n')
	print('Cannot find avconv executable, specify with "libavpath="')
	print('Cannot continue, aborting execution')
	exit(1)
if exiftool:
	if logtofile:
		logfile.write('Using exiftool found at %s\n\n' % which('exiftool'))
else:
	if logtofile:
		logfile.write('Cannnot find exiftool executable, falling back to defaults for folder names and duration\n\n')
	print('Cannot find exiftool executable, specify with "exiftoolpath="')
	print('Falling back to defaults for folder names and duration')
if args.moviepath == 'NULL':
	if logtofile:
		logfile.write('moviepath must not be NULL, exiting\n')
	print('moviepath must not be NULL. Stop')
	exit(1)
else:
	moviepath = args.moviepath
if args.imagepath == 'NULL':
	if logtofile:
		logfile.write('imagepath must not be NULL, exiting\n')
	print('imagepath must not be NULL. Stop')
	exit(1)
else:
	imagepath = args.imagepath
if args.x == 'NULL':
	if logtofile:
		logfile.write('x must not be NULL, exiting\n')
	print('x must not be NULL. Stop')
	exit(1)
else:
	x = args.x
if args.y == 'NULL':
	if logtofile:
		logfile.write('y must not be NULL, exiting\n')
	print('y must not be NULL. Stop')
	exit(1)
else:
	y = args.y
if args.fps == 'NULL':
	if logtofile:
		logfile.write('fps must not be NULL, exiting\n')
	print('fps must not be NULL. Stop')
	exit(1)
else:
	fps = args.fps
if args.nsec == 'NULL':
	nsec = None
else:
	nsec = float(args.nsec)
if args.start == 'NULL':
	start_given = None
else:
	start_given = float(args.start)
if args.stop == 'NULL':
	stop_given = None
else:
	stop_given = float(args.stop)
if args.ext == ['NULL']:
	if logtofile:
		logfile.write('ext must not be NULL, exiting\n')
	print('ext must not be NULL. Stop')
	exit(1)
else:
	extensions = []
	for e in args.ext:
		extensions.append(e.lower())
if args.verbose not in ['FALSE', 'TRUE']:
	if logtofile:
		logfile.write('Unknown value %s for verbose. Converting to FALSE\n\n')
	print('Unknown value %s for verbose. Converting to FALSE' % args.verbose)
	verbose = False
else:
	verbose = args.verbose == 'TRUE'
if start_given and stop_given and not start_given < stop_given:
	if logtofile:
		logfile.write('Error: start (%.2f) must be smaller than stop (%.2f), exiting\n' % (start_given, stop_given))
	print('Error: start (%.2f) must be smaller than stop (%.2f). Stop.' % (start_given, stop_given))
	exit(1)
if start_given and start_given < 0:
	if logtofile:
		logfile.write('Error: start must not be negative, exiting\n' % (start_given, stop_given))
	print('Error: start must not be negative. Stop.')
	exit(1)
if verbose:
	print('createImageSec %s' % version)
movieDir = os.path.abspath(moviepath)
sequenceParentDir = os.path.abspath(imagepath)	
if logtofile:
	logfile.write('Running script with values: localdir: %s, moviedir: ./%s, sequencedir: ./%s, x: %s, y: %s, fps: %s, nsec: %s, start_given: %s, stop_given: %s, extensions: %s\n\n' % (
		os.getcwd(), moviepath, imagepath, x, y, fps, nsec, start_given, stop_given, extensions
		)) 
def getDuration(filePath):
	if exiftool:
		try:
			durationRawFirst = subprocess.check_output([
				exiftool,
				'-api', 'LargeFileSupport',
				'-n',
				'-s3',
				'-duration',
				filePath
				], stderr=open('/dev/null')).split()[0]
			if python3:
#			  print('durationRawFirst: ', durationRawFirst, type(durationRawFirst))
			  duration = float(durationRawFirst.decode())
			else:
#			  print('durationRawFirst: ', durationRawFirst, type(durationRawFirst))
			  duration = float(str(durationRawFirst))
#			print(duration)
		except subprocess.CalledProcessError as e:
			if logtofile:
				logfile.write('Error in determining duration for file %s: %s\n\n' % (filePath, e))
			duration = None
	else:
		duration = None
	return duration
def getTargetDirNamePrefix(fileName):
	if logtofile:
		logfile.write('Entering getTargetDirNamePrefix with fileName = "%s"\n\n' % fileName)
	if exiftool:
		try:
			targetDirNamePrefixFromExif = subprocess.check_output([
				exiftool,
				'-api', 'LargeFileSupport',
				'-DateTimeOriginal',
				'-T',
				os.path.join(movieDir, fileName)
			], stderr=open('/dev/null'))
			if python3:
			  targetDirNamePrefixFromExif = targetDirNamePrefixFromExif.decode()
#			print(targetDirNamePrefixFromExif, type(targetDirNamePrefixFromExif))
			targetDirNamePrefixRaw = str(targetDirNamePrefixFromExif).strip().strip('-')
			# Note: '-' is sometimes returned by exiftool if DateTimeOriginal field empty

			if logtofile:
				logfile.write('Got targetDirNamePrefixRaw = "%s"\n\n' % targetDirNamePrefixRaw)
			if targetDirNamePrefixRaw == '':
				targetDirNamePrefix = nodate_prefix
			else:
				targetDirNamePrefixFirst = targetDirNamePrefixRaw.split()[0]
				targetDirNamePrefix = ''.join([c for c in targetDirNamePrefixFirst if c != ':']) + '_'
		except subprocess.CalledProcessError as e:
			if logtofile:
				logfile.write('Error in determining date for file %s: %s\n\n' % (os.path.join(movieDir, fileName), e))
				logfile.write('Got raw prefix: "%s"\n\n' % targetDirNamePrefixRaw)
			targetDirNamePrefix = nodate_prefix
		if not len(targetDirNamePrefix) == 8 and targetDirNamePrefix.isdigit():
			targetDirNamePrefix = nodate_prefix
		return targetDirNamePrefix
	else:
		return nodate_prefix
movieNames = []
if logtofile:
	logfile.write('Found the following entries in moviedir (%s):\n\n+++\n\n' % movieDir)
for fileName in os.listdir(movieDir):
	if logtofile:
		logfile.write('%s\n' % fileName)
	movieName, movieExtension = os.path.splitext(fileName)
	if not (os.path.isfile(os.path.join(movieDir, fileName))
			and movieExtension[1:].lower() in extensions):
		if logtofile:
			logfile.write('File %s has the wrong name or is a directory, skipping file\n\n' % fileName)
	else:
		movieNames.append(fileName)
if logtofile:
	logfile.write('\n+++\n\n')
if not os.path.exists(sequenceParentDir):
	os.mkdir(sequenceParentDir)
if start_given and stop_given and nsec:
	if logtofile:
		logfile.write('both start and stop given: nsec ignored\n\n')
	if verbose:
		print('Warning - both start and stop given: nsec ignored.')
for fileName in movieNames:
	movieName, movieExtension = os.path.splitext(fileName)
	targetDirNamePrefix = getTargetDirNamePrefix(fileName)
	targetDirName = targetDirNamePrefix + fileName
	if os.path.exists(os.path.join(sequenceParentDir, targetDirName)):
		if logtofile:
			logfile.write('Folder %s already exists, file %s skipped\n\n' % (targetDirName, fileName))		
		if verbose:
			print('Folder %s already exists, file %s skipped' % (targetDirName, fileName))
	else:
		if targetDirNamePrefix == nodate_prefix:
			if logtofile:
				logfile.write('File %s: exiftool cannot determine date\n\n' % fileName)		
			if verbose:
				print('File %s: exiftool cannot determine date' % fileName)
		if start_given:
			filestart = start_given
		else:
			filestart = 0
		if stop_given:
			filestop = stop_given
		else:
		  if python3:
		    filestop = sys.maxsize
		  else:
		    filestop = sys.maxint
		duration = getDuration(os.path.join(movieDir, fileName))
		if logtofile:
			logfile.write('Duration %s found for file %s\n\n' % (duration, fileName))
		if not duration:
			if logtofile:
				logfile.write('File %s: exiftool cannot determine length\n\n' % fileName)
			if verbose:
				print('File %s: exiftool cannot determine length.' % fileName)
		if start_given and stop_given:
			pass
		elif start_given:
			if nsec:
				filestop = filestart + nsec
		elif stop_given:
			if nsec:
				filestart = filestop - nsec
		else:
			if nsec:
				if duration:
				  filestart = duration/2 - nsec/2
				  filestop = duration/2 + nsec/2
				else:
					filestop = nsec
		if duration and not filestart < duration:
			if logtofile:
				logfile.write('File %s: Start (%.2f) should be smaller dan duration (%.2f sec). Skipping video\n\n' % (fileName, filestart, duration))
			if verbose:
				print('File %s: Start (%.2f) should be smaller dan duration (%.2f sec). Skipping video.' % (fileName, filestart, duration))
			continue
		filestart = max(filestart, 0)
		if duration and  filestop > duration:
			if logtofile:
				logfile.write('File %s: Requested end past end of file, setting end to %.2f\n\n' % (fileName, duration))
			if verbose:
				print('File %s: Requested end past end of file, setting end to %.2f' % (fileName, duration))
			filestop = min(filestop, duration)
		if filestop - filestart < nsec:
			if logtofile:
				logfile.write( 'File %s: Cannot convert requested duration (%.2f sec), converting %.2f sec\n\n' % (fileName, nsec, filestop-filestart))
			if verbose:
				print('File %s: Cannot convert requested duration (%.2f sec), converting %.2f sec.' % (fileName, nsec, filestop-filestart))
		command = conv_command(fileName, fps, x, y, targetDirName, filestart, filestop)
		if logtofile:
			logfile.write('Command: %s\n\n' % command)
		os.mkdir(os.path.join(sequenceParentDir, targetDirName))
		try:
			subprocess.call(command)
			if duration:
				if logtofile:
					logfile.write('File %s: Converting with start: %.2f, stop: %.2f to folder %s\n\n' % (fileName, filestart, filestop, targetDirName))
				if verbose:
					print('File %s: Converting with start: %.2f, stop: %.2f to folder %s.' % (fileName, filestart, filestop, targetDirName))
			else:
				if logfile:
					logfile.write('File %s: Attempting to convert with start: %.2f, stop: %.2f to folder %s\n\n' % (fileName, filestart, filestop, targetDirName))
				if verbose:
					print('File %s: Attempting to convert with start: %.2f, stop: %.2f to folder %s.' % (fileName, filestart, filestop, targetDirName))
		except Exception as e:
			logfile.write('Error in processing file %s: %s\n\n' % (fileName, e))
if logtofile:
	logfile.close()
