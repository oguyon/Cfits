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
function aoconflog {
echo "$@" >> aolconf.log
dolog "$LOOPNAME" "$@"
}

# external log, less verbose
function aoconflogext {
echo "$@" >> aolconf.log
dolog "$LOOPNAME" "$@"
dologext "$LOOPNAME $@"
}





function stringcenter {
line=$1
    let " col1 = $nbwcols-35"
    columns="$col1"
    string=$(printf "%*s%*s\n" $(( (${#line} + columns) / 2)) "$line" $(( (columns - ${#line}) / 2)) " ")
}






# DM DC level [um] for each Vmax setting

dmDCum025="0.0219"
dmDCum050="0.0875"
dmDCum075="0.1969"
dmDCum100="0.3500"
dmDMum125="0.5469"
dmDCum150="0.7857"




# Set DM Vmax 
function Set_dmVmax {
file="./status/status_dmVmax.txt"
currentfw=$(echo "$(cat $file)")
if [ ! "${current_dmVmax}" == "$1" ]
then
echo "CHANGING dmVmax to $1"# &> ${outmesg}

else
echo "WHEEL ALREADY SET TO"# > ${outmesg}
fi
sleep 0.1
currentfw=$1
echo "${currentfw}" > $file
sleep 0.1
}









state="menuhardwarecontrol"


while true; do

stateok=0

mkdir -p status


if [ $state = "menuhardwarecontrol" ]; then
stateok=1
menuname="HARDWARE CONTROL - LOOP ${LOOPNAME} ($LOOPNUMBER})\n"



file="./conf/conf_dmVmax.txt"
if [ -f $file ]; then
dmVmax=$(cat $file)
else
dmVmax="125"
echo "$dmVmax" > $file
fi

file="./conf/conf_dmDCum.txt"
if [ -f $file ]; then
dmDCum=$(cat $file)
else
dmDCum="125"
echo "$dmDCum" > $file
fi



stringcenter "DM control  [ current: Vmax = ${dmVmax} V  DC = ${dmDCum} um ]"
menuitems=( "1 ->" "\Zb\Zr$string\Zn" )
menuitems+=( "" "" )

menuitems+=( "" "" )


if [ "$dmVmax" = " 25" ]; then
menuitems+=( "dmVmax025" "\Zr\Z2 dmVmax =  25 V  (DC level = ${dmDCum025} um)\Zn" )
else
menuitems+=( "dmVmax025" " dmVmax =  25 V  (DC level = ${dmDCum025} um)" )
fi

if [ "$dmVmax" = " 50" ]; then
menuitems+=( "dmVmax050" "\Zr\Z2 dmVmax =  50 V  (DC level = ${dmDCum050} um)\Zn" )
else
menuitems+=( "dmVmax050" " dmVmax =  50 V  (DC level = ${dmDCum050} um)" )
fi

if [ "$dmVmax" = " 75" ]; then
menuitems+=( "dmVmax075" "\Zr\Z2 dmVmax =  75 V  (DC level = ${dmDCum075} um)\Zn" )
else
menuitems+=( "dmVmax075" " dmVmax =  75 V  (DC level = ${dmDCum075} um)" )
fi

if [ "$dmVmax" = "100" ]; then
menuitems+=( "dmVmax100" "\Zr\Z2 dmVmax = 100 V  (DC level = ${dmDCum100} um)\Zn" )
else
menuitems+=( "dmVmax100" " dmVmax = 100 V  (DC level = ${dmDCum100} um)" )
fi

if [ "$dmVmax" = "125" ]; then
menuitems+=( "dmVmax125" "\Zr\Z2 dmVmax = 125 V  (DC level = ${dmDCum125} um)\Zn" )
else
menuitems+=( "dmVmax125" " dmVmax = 125 V  (DC level = ${dmDCum125} um)" )
fi

if [ "$dmVmax" = "150" ]; then
menuitems+=( "dmVmax150" "\Zr\Z2 dmVmax = 150 V  (DC level = ${dmDCum150} um)\Zn" )
else
menuitems+=( "dmVmax150" " dmVmax = 150 V  (DC level = ${dmDCum150} um)" )
fi








stringcenter "Science IRcam"
menuitems+=( "2 ->" "\Zb\Zr$string\Zn" )
menuitems+=( "" "" )




dialog --colors --title "Alignment" \
--ok-label "Select" \
--cancel-label "Top" \
--help-button --help-label "Exit" \
--default-item "${menuhardwarecontrol_default}" \
 --menu "$menuname" \
  $nbwlines $nbwcols 100 "${menuitems[@]}"  2> $tempfile


retval=$?
choiceval=$(cat $tempfile)


menualign_default="$choiceval"
state="menuhardwarecontrol"




case $retval in
   0) # button
menucontrolhardware_default="$choiceval"
	case $choiceval in


	dmVmax025)
dmVmax=" 25"
echo "${dmVmax}" > ./conf/conf_dmVmax.txt
Cfits << EOF
aolsetdmvoltmax 00 ${dmVmax}
aolsetdmDC 00 ${dmDCum025}
exit
EOF
;;

	dmVmax050)
dmVmax=" 50"
echo "${dmVmax}" > ./conf/conf_dmVmax.txt
Cfits << EOF
aolsetdmvoltmax 00 ${dmVmax}
aolsetdmDC 00 ${dmDCum050}
exit
EOF
;;

	dmVmax075)
dmVmax=" 75"
echo "${dmVmax}" > ./conf/conf_dmVmax.txt
Cfits << EOF
aolsetdmvoltmax 00 ${dmVmax}
aolsetdmDC 00 ${dmDCum075}
exit
EOF
;;

	dmVmax100)
dmVmax="100"
echo "${dmVmax}" > ./conf/conf_dmVmax.txt
Cfits << EOF
aolsetdmvoltmax 00 ${dmVmax}
aolsetdmDC 00 ${dmDCum100}
exit
EOF
;;

	dmVmax125)
dmVmax="125"
echo "${dmVmax}" > ./conf/conf_dmVmax.txt
Cfits << EOF
aolsetdmvoltmax 00 ${dmVmax}
aolsetdmDC 00 ${dmDCum125}
exit
EOF
;;

	dmVmax150)
dmVmax="150"
echo "${dmVmax}" > ./conf/conf_dmVmax.txt
Cfits << EOF
aolsetdmvoltmax 00 ${dmVmax}
aolsetdmDC 00 ${dmDCum150}
exit
EOF
;;




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
