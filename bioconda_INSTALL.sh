#!/usr/bin/env bash

set -e # Exit as soon as any line in the bash script fails

ROOTDIR=$( cd $(dirname $0) ; pwd -P ) # path to main PiReT directory

echo
exec &> >(tee -a  install.log)
exec 2>&1 # copies stderr onto stdout

# create a directory where all dependencies will be installed
cd $ROOTDIR
mkdir -p thirdParty
cd thirdParty


utility_tools=(samtools bedtools hisat2)
alignments_tools=(bowtie2 bwa)
perl_modules=( perl_parallel_forkmanager )
all_tools=("${utility_tools[@]}" "${alignments_tools[@]}")


install_hisat2()
{
echo "--------------------------------------------------------------------------
                           installing hisat2 v2.0.5
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda hisat2=2.0.5
echo "
------------------------------------------------------------------------------
                           hisat2 v2.0.5 installed
------------------------------------------------------------------------------
"
}

install_jellyfish()
{
echo "--------------------------------------------------------------------------
                           installing jellyfish v2.2.6
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda jellyfish=2.2.6
echo "
------------------------------------------------------------------------------
                           jellyfish v2.2.6 installed
------------------------------------------------------------------------------
"
}

install_perl_parallel_forkmanager()
{
echo "------------------------------------------------------------------------------
               Installing Perl Module Parallel-ForkManager v1.17
------------------------------------------------------------------------------
"
conda install --yes -c bioconda perl-parallel-forkmanager=1.17

echo "
------------------------------------------------------------------------------
                        Parallel-ForkManager-1.17 Installed
------------------------------------------------------------------------------
"
}

install_bowtie2()
{
echo "--------------------------------------------------------------------------
                           installing bowtie2 v2.2.8
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda bowtie2=2.2.8
echo "
------------------------------------------------------------------------------
                           bowtie2 v2.2.8 installed
------------------------------------------------------------------------------
"
}

install_bwa()
{
echo "------------------------------------------------------------------------------
                           Downloading bwa v0.7.15
------------------------------------------------------------------------------
"
conda install --yes -c bioconda bwa=0.7.15
echo "
------------------------------------------------------------------------------
                           bwa v0.7.15 installed
------------------------------------------------------------------------------
"
}

install_htseq()
{
echo "------------------------------------------------------------------------------
                           downloading htseq v0.6.1p1
------------------------------------------------------------------------------
"
conda install --yes -c bioconda htseq
echo "
------------------------------------------------------------------------------
                           htseq v0.6.1p1 installed
------------------------------------------------------------------------------
"
}


install_samtools()
{
echo "--------------------------------------------------------------------------
                           Downloading samtools v1.3.1
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda samtools=1.3.1
echo "
--------------------------------------------------------------------------------
                           samtools v1.3.1 installed
--------------------------------------------------------------------------------
"
}

install_cpan()
{
echo "--------------------------------------------------------------------------
                           Installing cpan
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda perl-cpan-meta
echo "
--------------------------------------------------------------------------------
                           cpan installed
--------------------------------------------------------------------------------
"
}


install_bedtools()
{
echo "--------------------------------------------------------------------------
                           Compiling bedtools v2.26.0
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda bedtools=2.26.0
echo "
--------------------------------------------------------------------------------
                           bedtools v2.26.0 compiled
--------------------------------------------------------------------------------
"
}

install_miniconda()
{
echo "--------------------------------------------------------------------------
                           downloading miniconda
--------------------------------------------------------------------------------
"

if [[ "$OSTYPE" == "darwin"* ]]
then
{

  curl -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda2-4.2.12-MacOSX-x86_64.sh
  chmod +x miniconda.sh
  ./miniconda.sh -b -p $ROOTDIR/thirdParty/miniconda -f
  export PATH=$ROOTDIR/thirdParty/miniconda/bin:$PATH

}
else
{  

  wget https://repo.continuum.io/miniconda/Miniconda2-4.2.12-Linux-x86_64.sh -O miniconda.sh
  chmod +x miniconda.sh
  ./miniconda.sh -b -p $ROOTDIR/thirdParty/miniconda -f
  export PATH=$ROOTDIR/thirdParty/miniconda/bin:$PATH

}
fi
}

install_perl_string_approx()
{
echo "--------------------------------------------------------------------------
                           installing Perl Module String::Approx
--------------------------------------------------------------------------------
"
cpan String::Approx
echo "
--------------------------------------------------------------------------------
                           String::Approx installed
--------------------------------------------------------------------------------
"
}

install_gffread()
{
echo "--------------------------------------------------------------------------
                      installing gffread
--------------------------------------------------------------------------------
"
cd $ROOTDIR/thirdParty/miniconda/bin
git clone https://github.com/gpertea/gclib
git clone https://github.com/gpertea/gffread gffread_git
cd gffread_git
make
cp gffread ../
echo "
--------------------------------------------------------------------------------
                           gffread installed
--------------------------------------------------------------------------------
"
}



checkSystemInstallation()
{
    IFS=:
    for d in $PATH; do
      if test -x "$d/$1"; then return 0; fi
    done
    return 1
}

checkLocalInstallation()
{
    IFS=:
    for d in $ROOTDIR/bin; do
      if test -x "$d/$1"; then return 0; fi
    done
    return 1
}

checkPerlModule()
{
   perl -e "use lib \"$rootdir/lib\"; use $1;"
   return $?
}


containsElement () {
  local e
  for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
  return 1
}

print_usage()
{
cat << EOF
usage: $0 options
    If no options, it will check existing installation and run tools installation for those uninstalled.
    options:
    help            show this help
    list            show available tools for updates
    tools_name      install/update individual tool
    force           force to install all list tools locally
    
    ex: To update bowtie2 only
        $0 bowtie2
    ex: To update bowtie2 and bwa
        $0 bowtie2 bwa
    ex: RE-install Phylogeny tools
        $0 Phylogeny
        
EOF

}

print_tools_list()
{


   echo "Available tools for updates/re-install"
   echo -e "\nAlignment"
   for i in "${alignments_tools[@]}"
   do
       echo "* $i"
   done
   echo -e "\nUtility"
   for i in "${utility_tools[@]}"
   do
       echo "* $i"
   done
}


### Main ####
if ( checkSystemInstallation csh )
then
  #echo "csh is found"
  echo ""
else
  echo "csh is not found"
  echo "Please Install csh first, then INSTALL the package"
  exit 1
fi

if [ "$#" -ge 1 ]
then
  for f in $@
  do
    case $f in
      help)
        print_usage
        exit 0;;
      list)
        print_tools_list
        exit 0;;
      Alignment)
        for tool in "${alignments_tools[@]}"
        do
            install_$tool
        done
        echo -e "Alignment tools installed.\n"
        exit 0;;
      Utility)
        for tool in "${utility_tools[@]}"
        do
            install_$tool
        done
        echo -e "Utility tools installed.\n"
        exit 0;;
      force)
        for tool in "${all_tools[@]}"
        do
            install_$tool
        done
        ;;
      *)
        if ( containsElement "$f" "${alignments_tools[@]}" || containsElement "$f" "${utility_tools[@]}" )
        then
            install_$f
        else
            echo "$f: no this tool in the list"
            print_tools_list
        fi
        exit 0;;
    esac
  done
fi


if [[ "$OSTYPE" == "darwin"* ]]
then
{
    if ( checkSystemInstallation R )
    then
    {
        echo "R is found"
    }
    else
    {
        echo "R is not found"
        echo "Please install R from http://cran.r-project.org/bin/macosx/";
        exit 1
    }
    fi
}
else
{  #TODO: Ask chien-chi why local installation is checked here
    if ( checkSystemInstallation R )
    then
    {
        echo "R is found"
    }
    else
    {
        echo "R is not found"
        install_R
    }
    fi
}
fi

# check if required bioconductor R packages are installed

#edgeR
echo "if(\"edgeR\" %in% rownames(installed.packages()) == FALSE)  {source('https://bioconductor.org/biocLite.R')
      biocLite('edgeR')}" | Rscript -

#Deseq2
echo "if(\"edgeR\" %in% rownames(installed.packages()) == FALSE)  {source('https://bioconductor.org/biocLite.R')
      biocLite('DESeq2')}" | Rscript -


if ( checkSystemInstallation conda )
then
  echo "conda is found"
else
  echo "conda is not found"
  install_miniconda
fi

if ( checkSystemInstallation hisat2 )
then
  echo "hisat2 is found"
else
  echo "hisat2 is not found"
  install_hisat2
fi

if ( checkSystemInstallation htseq-count )
then
  echo "htseq is found"
else
  echo "htseq is not found"
  install_htseq
fi

if ( checkSystemInstallation jellyfish )
then
  echo "jellyfish is found"
else
  echo "jellyfish is not found"
  install_jellyfish
fi

if ( checkSystemInstallation bowtie2 )
then
  echo "bowtie2 is found"
else
  echo "bowtie2 is not found"
  install_bowtie2
fi

if ( checkSystemInstallation bwa )
then
  echo "bwa is found"
else
  echo "bwa is not found"
  install_bwa
fi

if ( checkSystemInstallation samtools )
then
  echo "samtools is found"
else
  echo "samtools is not found"
  install_samtools
fi

if ( checkPerlModule Parallel::ForkManager )
then
  echo "Perl Parallel::ForkManager is found"
else
  echo "Perl Parallel::ForkManager is not found"
  install_perl_parallel_forkmanager
fi

if ( checkSystemInstallation bedtools )
then
  echo "bedtools is found"
else
  echo "bedtools is not found"
  install_bedtools
fi

if ( checkSystemInstallation cpan )
then
  echo "cpan is found"
else
  echo "cpan is not found"
  install_cpan
fi

if ( checkSystemInstallation gffread )
then
  echo "gffread is found"
else
  echo "gffread is not found"
  install_gffread
fi

if ( checkPerlModule String::Approx )
then
  echo "Perl String::Approx is found"
else
  echo "Perl String::Approx is not found"
  install_perl_string_approx
fi


################################################################################
#                       Add path to bash
################################################################################
if [ -f $HOME/.bashrc ]
then
{
  echo "#Added by RNASeq pipeline installation" >> $HOME/.bashrc
  echo "export RNASeq_HOME=$ROOTDIR" >> $HOME/.bashrc
  echo "export PATH=$ROOTDIR/thirdParty/miniconda/bin/:$PATH:$ROOTDIR/scripts" >> $HOME/.bashrc
  source $HOME/.bashrc 
}
else
{
  echo "#Added by RNASeq pipeline installation" >> $HOME/.bash_profile
  echo "export RNASeq_HOME=$ROOTDIR" >> $HOME/.bash_profile
  echo "export PATH=$ROOTDIR/thirdParty/miniconda/bin/:$PATH:$ROOTDIR/scripts" >> $HOME/.bash_profile
  source $HOME/.bash_profile 
}
fi

echo "
All done! Please Restart the Terminal Session.
Run
./runPipeline_rRNA_noqsub_commandline.pl
for usage.
Read the README for more information!
Thanks!
    "
