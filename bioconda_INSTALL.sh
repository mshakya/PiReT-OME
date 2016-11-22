#!/usr/bin/env bash
set -e
rootdir=$( cd $(dirname $0) ; pwd -P )
exec &> >(tee -a  install.log)
exec 2>&1

cd $rootdir




utility_tools=(samtools bedtools)
alignments_tools=( bowtie2 bwa )
all_tools=("${utility_tools[@]}" "${alignments_tools[@]}")


install_bowtie2()
{
echo "--------------------------------------------------------------------------
                           Compiling bowtie2
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda bowtie2
echo "
------------------------------------------------------------------------------
                           bowtie2 installed
------------------------------------------------------------------------------
"
}

install_bwa()
{
echo "------------------------------------------------------------------------------
                           Compiling bwa 
------------------------------------------------------------------------------
"
conda install --yes -c bioconda bwa
echo "
------------------------------------------------------------------------------
                           bwa installed
------------------------------------------------------------------------------
"
}


install_samtools()
{
echo "--------------------------------------------------------------------------
                           Compiling samtools 0.1.19
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda samtools
echo "
--------------------------------------------------------------------------------
                           samtools compiled
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
    for d in $rootdir/bin; do
      if test -x "$d/$1"; then return 0; fi
    done
    return 1
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

# echo "if(\"gridExtra\" %in% rownames(installed.packages()) == FALSE)  {install.packages(\"gridExtra_0.9.1.tar.gz\", repos = NULL, type=\"source\")}" | Rscript -  

#edgeR
echo "if(\"edgeR\" %in% rownames(installed.packages()) == FALSE)  {source('https://bioconductor.org/biocLite.R')
      biocLite('edgeR')}" | Rscript -  



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


if [ -f $HOME/.bashrc ]
then
{
  echo "#Added by RNASeq pipeline installation" >> $HOME/.bashrc
  echo "export RNASeq_HOME=$rootdir" >> $HOME/.bashrc
  echo "export PATH=$rootdir/bin/:$PATH:$rootdir/scripts" >> $HOME/.bashrc
}
else
{
  echo "#Added by RNASeq pipeline installation" >> $HOME/.bash_profile
  echo "export RNASeq_HOME=$rootdir" >> $HOME/.bash_profile
  echo "export PATH=$rootdir/bin/:$PATH:$rootdir/scripts" >> $HOME/.bash_profile
}
fi

# sed -i.bak 's,%EDGE_HOME%,'"$rootdir"',g' $rootdir/edge_ui/cgi-bin/edge_config.tmpl
# sed -i.bak 's,%EDGE_HOME%,'"$rootdir"',g' $rootdir/edge_ui/apache_conf/edge_apache.conf

# TOLCPU=`cat /proc/cpuinfo | grep processor | wc -l`;
# if [ $TOLCPU -gt 0 ]
# then
# {
#     sed -i.bak 's,%TOTAL_NUM_CPU%,'"$TOLCPU"',g' $rootdir/edge_ui/cgi-bin/edge_config.tmpl
#     DEFAULT_CPU=`echo -n $((TOLCPU/3))`;
#     if [ $DEFAULT_CPU -lt 1 ]
#     then
#     {
#         sed -i.bak 's,%DEFAULT_CPU%,'"1"',g' $rootdir/edge_ui/index.html
#     }
#     else
#     {
#         sed -i.bak 's,%DEFAULT_CPU%,'"$DEFAULT_CPU"',g' $rootdir/edge_ui/index.html
#     }
#     fi
# }
# fi

echo "
All done! Please Restart the Terminal Session.
Run
./runPipeline
for usage.
Read the README for more information!
Thanks!
"
