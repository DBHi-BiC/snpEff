

									Release process
									---------------

Run JUnit tests (TestSuiteAll) and make sure all of them pass

Run integration tests

Create JAR files

	cd ~/snpEff
	./scripts/make.sh 

Download ENSEMBL databases:

	./scripts/download.sh

	# Copy queue script to databases
	cp queue_build.txt scripts/queue_build.txt

Build databases:

	# Edit script to adapt for max number of processes 
	./scripts/queue_build.sh

Check that databases have been built correctly:

	grep Error *.stdout | cut -f 2- -d : | sort > check.txt
	
	cat check.txt | grep -v "0.0%" > check.non_zero.txt 
	vi check.non_zero.txt
	# Note: All error percentages should be below 1%

- Update download page

- Update Galaxy's snpEff.xml

- Create a package snpEff_VVV.zip 

	./scripts/distro.sh

- Upload to sourceforge
		
		# Core program
		scp snpEff_v3_1_core.zip pcingola,snpeff@frs.sourceforge.net:/home/frs/project/s/sn/snpeff/snpEff_latest_core.zip
		scp snpEff_v3_1_core.zip pcingola,snpeff@frs.sourceforge.net:/home/frs/project/s/sn/snpeff/
				
		# Individual databases
		scp snpEff_v3_1_*.zip pcingola,snpeff@frs.sourceforge.net:/home/frs/project/s/sn/snpeff/databases/v3_1/
		
		# SnpSift
		scp SnpSift.jar pcingola,snpeff@frs.sourceforge.net:/home/frs/project/s/sn/snpeff/SnpSift_v1_7.jar
		scp SnpSift.jar pcingola,snpeff@frs.sourceforge.net:/home/frs/project/s/sn/snpeff/SnpSift_latest.jar
		
		
- Update snpEff pages (at least download link). 
		
		# Create download table (add output to download_content.html)
		java -jar snpEff.jar cfg2table download 

		# Create galaxy menu (add output to galaxy/snpEff.xml)
		java -jar snpEff.jar cfg2table galaxy 
		
		cd $HOME/workspace/SnpEff/html/
		
		# Copy html files 
		scp style.css *.html pcingola,snpeff@frs.sourceforge.net:htdocs/
		
		# Copy images
		scp -r  images/ pcingola,snpeff@frs.sourceforge.net:htdocs/images/

- Upload to Galaxy ToolShed: http://toolshed.g2.bx.psu.edu/
		Reference: http://wiki.g2.bx.psu.edu/Tool%20Shed
		
		
-------------------------------------------------------------------------------		
		
How to add the libraries to your local Maven repository:

    http://maven.apache.org/general.html#importing-jars

    mvn install:install-file -Dfile=picard-1.77.jar -DgroupId=net.sf.picard -DartifactId=Picard -Dversion=1.77 -Dpackaging=jar
    mvn install:install-file -Dfile=sam-1.77.jar -DgroupId=net.sf.samtools -DartifactId=Sam -Dversion=1.77 -Dpackaging=jar
   
