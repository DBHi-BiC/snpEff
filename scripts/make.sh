#!/bin/sh

VERSION_SNPEFF=3.1


#---
# Build SnpEff
#---

cd $HOME/workspace/SnpEff/

mvn assembly:assembly
cp target/snpEff-$VERSION_SNPEFF-jar-with-dependencies.jar $HOME/workspace/SnpEff/snpEff.jar


# Install JAR file in local Maven repo
mvn install:install-file \
	-Dfile=target/snpEff-$VERSION_SNPEFF.jar \
	-DgroupId=ca.mcgill.mcb.pcingola \
	-DartifactId=snpEff \
	-Dversion=$VERSION_SNPEFF \
	-Dpackaging=jar \
	-DgeneratePom=true

cd - 
