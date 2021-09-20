#!/bin/bash

logfile="log.txt" #creates logfile called log.txt when used in log() function 

log(){
        echo "`date +"%F %T"` ${1}" >> ${logfile} #passes date, time, and messages to log.txt
        echo "`date +"%F %T"` ${1}" 
}

parameter_check(){
        if [ -z $1 ]; then #if there is not a first parameter passed to code then
                log "Parameter 1 missing" # print message to log.txt
                exit 1 #missing parameter returns exit status of 1
        fi
}

create_folder(){
        if [ ! -d $1 ]; then #if there is not a directory named as the parameter passed to create_folder then
                log "making dir: $1" #print message to log.txt
                mkdir $1 #creare directory named after parameter 1 
        fi
}

folder_check(){
        if [ -d $1 ]; then #if there is a directory named after the parameter passed to folder_check then
        log "directory $1 already exists" #print message to log.txt
        exit
        fi
        }


file_check(){
        if [ -e $1 ]; then #if a file named after parameter one is found then 
        log "$1 already exists" #print message to log.txt
        exit
        fi
}
