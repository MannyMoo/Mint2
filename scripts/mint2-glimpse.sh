#!/bin/bash

find $MINT2ROOT/src $MINT2ROOT/Mint -type f -exec grep -H "$1" {} \;
