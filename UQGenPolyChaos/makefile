# type 'make -n' to see what this will do
MovieFiles = linearperfusiondog1.gif \
             linearconductivitydog1.gif \
             nonlinearperfusiondog1.gif \
             linearabsorptiondog1.gif \
             linearscatteringdog1.gif \
             linearscatteringabsorptiondog1.gif \
             linearscatteringabsorptionperfusiondog1.gif \
             linearscatteringabsorptionnonlinearperfusiondog1.gif

# type 'make --just-print' to see what this will do
all: $(MovieFiles)
	for file in $(MovieFiles); do \
	$(MAKE) $$file; \
	done

#default values
DELAY = 30
FILESKIP = 1 # dont skip any files

# convert avi to gif
# "FORCE" doesn't seem to be affecting automatic pattern rules
#         http://www.gnu.org/software/make/manual/make.html#Pattern-Rules
#echo -e  "import os\nfor i in range(1,len(os.listdir('$*'))):\n  if(i%$(FILESKIP) > 0):   os.remove('$*/%08d.jpg' % i) " | python
%.gif: %.avi FORCE
	mplayer $< -vo jpeg:outdir=$*
	convert -delay $(DELAY) $(RESIZE) -loop 0 $*/000000[45678]*.jpg $@

#http://www.gnu.org/software/make/manual/make.html#Force-Targets
FORCE:


