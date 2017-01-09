#!/bin/bash

execname="Cfits"


tempfile=`tempfile 2>/dev/null` || tempfile=/tmp/test$$
trap "rm -f $tempfile" 0 1 2 5 15


LINES=$( tput lines )
COLUMNS=$( tput cols )
let " nbwlines = $LINES - 10 "
let " nbwcols = $COLUMNS - 10 "
echo "$COLUMNS -> $nbwcols"





LOOPNUMBER_file="LOOPNUMBER"
confnbfile="./conf/conf_CONFNUMBER.txt"


mkdir -p conf
mkdir -p status
mkdir -p tmp


LOOPNUMBER_default=2  # loop number

# LOOPNUMBER (loop number)
if [ ! -f $LOOPNUMBER_file ]
then
	echo "creating loop number"
	echo "$LOOPNUMBER_default" > $LOOPNUMBER_file
else
	LOOPNUMBER=$(cat $LOOPNUMBER_file)
	echo "LOOPNUMBER = $LOOPNUMBER"
fi




# ======================= LOGGING =================================
LOOPNAME=$( cat LOOPNAME )
echo "LOOPNAME = $LOOPNAME"
# internal log - logs EVERYTHING
function FPaoconflog {
echo "$@" >> FPaolconf.log
dolog "$LOOPNAME" "$@"
}

# external log, less verbose
function FPaoconflogext {
echo "$@" >> FPaolconf.log
dolog "$LOOPNAME" "$@"
dologext "$LOOPNAME $@"
}





function stringcenter {
line=$1
    let " col1 = $nbwcols-35"
    columns="$col1"
    string=$(printf "%*s%*s\n" $(( (${#line} + columns) / 2)) "$line" $(( (columns - ${#line}) / 2)) " ")
}









state="menualign"


while true; do

stateok=0

mkdir -p status


statfile="./status/status_alignPcam.txt"
Pcamloopstat=$(cat $statfile)
if [[ -f "$statfile" && ( "$Pcamloopstat" = " ON" || "$Pcamloopstat" = "OFF" || "$Pcamloopstat" = "PAU" ) ]]; then
echo "OK"
else
echo "OFF" > $statfile
Pcamloopstat="OFF"
fi


PyrFilter=$(cat ./status/status_fw.txt)


if [ $state = "menualign" ]; then
stateok=1
menuname="ALIGNMENT - LOOP ${LOOPNAME} ($LOOPNUMBER})\n
\n"





stringcenter "Alignment"
menuitems=( "1 ->" "\Zb\Zr$string\Zn" )




menuitems+=( "" "" )


menuitems+=( "" "" )







dialog --colors --title "Alignment" \
--ok-label "Select" \
--cancel-label "Top" \
--help-button --help-label "Exit" \
--default-item "${menualign_default}" \
 --menu "$menuname" \
  $nbwlines $nbwcols 100 "${menuitems[@]}"  2> $tempfile


retval=$?
choiceval=$(cat $tempfile)


menualign_default="$choiceval"
state="menualign"

case $retval in
   0) # button
menualign_default="$choiceval"
	case $choiceval in



	esac;;
   1) state="menutop";;   
   2) state="menuexit";;
   255) state="menuexit";;
esac
fi













if [ $state = "menuexit" ]; then
stateok=1
echo "menuexit -> exit"
exit
fi



if [ $stateok = 0 ]; then
echo "state \"$state\" not recognized ... exit"
exit
fi




done
