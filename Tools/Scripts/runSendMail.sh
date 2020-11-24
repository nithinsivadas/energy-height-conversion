#!/bin/bash

export MAILTO=$1
export SUBJECT=$2
export BODY=$3
export ATTACH=$4
export SENDER="nithin@bu.edu"
export NAME="Nithin Sivadas"
(
 echo "To: $MAILTO"
 echo "Subject: $SUBJECT"
 echo "From: $NAME <$SENDER>"
 echo "MIME-Version: 1.0"
 echo 'Content-Type: multipart/mixed; boundary="-q1w2e3r4t5"'
 echo
 echo '---q1w2e3r4t5'
 echo "Content-Type: text/plain"
 echo "Content-Disposition: inline"
 cat $BODY
 echo '---q1w2e3r4t5'
 echo 'Content-Type: application; name="'$(basename $ATTACH)'"'
 echo "Content-Transfer-Encoding: base64"
 echo 'Content-Disposition: attachment; filename="'$(basename $ATTACH)'"'
 uuencode --base64 $ATTACH $(basename $ATTACH)
 echo '---q1w2e3r4t5--'
) | sendmail $MAILTO
