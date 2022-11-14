#set -x
#SCRIPT PARA AUXILIAR NA EXECUCAO DE PROGRAMAS 
# DE SIMULACAO DE ESCOAMENTOS EM RESERVATORIOS 

#### EXPERIMENTOS ####

# diretorio do experimento:
# dirPadrao ou lido de primeiro argumento da linha de comando
dirPadrao="experiment"
dirExp="$(pwd)"/${1:-${dirPadrao}} 
arqTela="$(pwd)"/${2:-${dirExp}/tela.txt} 

rm -rf ${dirExp}/out/*.vtk
rm -rf ${dirExp}/out/*.vtu
rm -rf ${dirExp}/fort.*

#### DEFINICAO DO NUMERO DE THREADS (OPENMP) E PROCESSOS (MPI) ####
ntPadrao=1
numThreads=${3:-${ntPadrao}}
export OMP_NUM_THREADS=${numThreads}


npPadrao=1
#numProcs=${4:-${npPadrao}}
numProcs=${npPadrao}
export NP=${numProcs} # atribuir valor padrao aa variavel NP

#ARGSP=" "
#echo "${ARGSP}"

#### DEFINICAO DO EXECUTAVEL ####
LOCAL=$(pwd)
NOMEEXECUTAVEL=simulador.exe
DIRBIN="${LOCAL}/bin"

EXECUTAVEL=${DIRBIN}/${NOMEEXECUTAVEL}
arqTela="./out/tela.txt"

#### definicao do comando a ser executado
# comando="(export OMP_NUM_THREADS=${numThreads} ; cd ${dirExp}; time ${EXECUTAVEL} |tee -a ${arqTela})"
comando="(export OMP_NUM_THREADS=${numThreads} ; cd ${dirExp}; time ${EXECUTAVEL} |tee ${arqTela})"

if [ -e ${EXECUTAVEL} ] 
then
  printf "\n diretorio do experimento.: %s\n" ${dirExp}  
  printf "\n nome do executavel.......: %s\n" ${EXECUTAVEL} 
  printf "\n numero de threads .......: %d\n" ${OMP_NUM_THREADS}
  printf "\n numero de processos......: %d\n" ${NP}
  printf "\n comando .................: %s\n" "${comando}"
  eval ${comando} # |tee  ${arqTela}
else
  printf "\n EXECUTAVEL NAO ENCONTRADO \n"
  printf "\n comando .................: %s\n" "${comando}"
fi
