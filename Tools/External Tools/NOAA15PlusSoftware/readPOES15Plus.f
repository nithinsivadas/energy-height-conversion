      program readPOES15Plus

      ! The program takes two arguments, specifying the input
      ! filename and the output filename (including paths.)
      ! The output file is written in the current working directory
      ! unless a path is included as part of the filename.

      integer*2 j,k,l
      integer*4 rec,eof,hr,min,sec
      character*80 arcfile,outfile
      ! Common block /rec/ variables:
      integer*4 cSum,ihd(6,4)
      integer*2 cSumFlag,major,status(10),qual(16),minor(16),mdf(40,16)
      real*4 analog(17),head(27,4),ssLoc(2,16),mep0(9,16),mep90(9,16),
     +    mepOmni(4,16),ted0(8,16),ted30(8,16),ted0s(8,4),ted30s(8,4),
     +    tedback(2,4),tedfx(7,16)

      common /rec/cSum,ihd,cSumFlag,major,status,qual,minor,
     +    mdf,analog,head,ssLoc,mep0,mep90,mepOmni,ted0,ted30,
     +    ted0s,ted30s,tedback,tedfx

      eof = 0
      call getarg(1,arcfile)
      call getarg(2,outfile)
      open(unit=7,access='direct',recl=2544,file = arcfile,status='old')
      open(unit=11,access='sequential',file = outfile,
     +     status='unknown',form='formatted')

      write(11,300)"year doy hr min sec lat lon alt(km) orbit "//
     +   "mep0P1 mep0P2 mep0P3 mep0P4 mep0P5 mep0P6 mep0E1 mep0E2 "//
     +   "mep0E3 mep90P1 mep90P2 mep90P3 mep90P4 "//
     +   "mep90P5 mep90P6 mep90E1 mep90E2 mep90E3 "//
     +   "mepOmniP6 mepOmniP7 mepOmniP8 mepOmniP9 "//
     +   "ted01 ted02 ted03 ted04 ted05 ted06 ted07 ted08 "//
     +   "ted301 ted302 ted303 ted304 ted305 ted306 ted307 ted308 "//
     +   "ted0s1 ted0s2 ted0s3 ted0s4 "//
     +   "ted0s5 ted0s6 ted0s7 ted0s8 "//
     +   "ted30s1 ted30s2 ted30s3 ted30s4 "//
     +   "ted30s5 ted30s6 ted30s7 ted30s8 "//
     +   "tedback11 tedback12 tedback13 tedback14 "//
     +   "tedback21 tedback22 tedback23 tedback24 "//
     +   "tedfx1 tedfx2 tedfx3 tedfx4 tedfx5 tedfx6 tedfx7"

      do rec = 1,2700
          call readArcSub(7,rec,eof)
          if(eof.eq.1) goto 100
          l = 1
          do j=1,4
              hr = ihd(4,j)/3600000
              min = (ihd(4,j) - hr*3600000)/60000
              sec = (ihd(4,j) - hr*3600000 - min*60000)/1000
              do k=1,4
              write(11,200),ihd(2,j),ihd(3,j),hr,min,sec,
     +                  head(2,j),head(3,j),ihd(5,j),ihd(6,j),
     +                  mep0(1,l),mep0(2,l),mep0(3,l),mep0(4,l),
     +          mep0(5,l),mep0(6,l),mep0(7,l),mep0(8,l),mep0(9,l),
     +          mep90(1,l),mep90(2,l),mep90(3,l),
     +          mep90(4,l),mep90(5,l),mep90(6,l),mep90(7,l),mep90(8,l),
     +          mep90(9,l),
     +          mepOmni(1,l),mepOmni(2,l),mepOmni(3,l),mepOmni(4,l),
     +          ted0(1,l),ted0(2,l),ted0(3,l),ted0(4,l),
     +          ted0(5,l),ted0(6,l),ted0(7,l),ted0(8,l),
     +          ted30(1,l),ted30(2,l),ted30(3,l),ted30(4,l),
     +          ted30(5,l),ted30(6,l),ted30(7,l),ted30(8,l),
     +          ted0s(1,j),ted0s(2,j),ted0s(3,j),
     +          ted0s(4,j),ted0s(5,j),ted0s(6,j),ted0s(7,j),ted0s(8,j),
     +          ted30s(1,j),ted30s(2,j),ted30s(3,j),ted30s(4,j),
     +          ted30s(5,j),ted30s(6,j),ted30s(7,j),ted30s(8,j),
     +          tedback(1,1),tedback(1,2),tedback(1,3),tedback(1,4),
     +          tedback(2,1),tedback(2,2),tedback(2,3),tedback(2,4),
     +          tedfx(1,j),tedfx(2,j),tedfx(3,j),tedfx(4,j),tedfx(5,j),
     +          tedfx(6,j),tedfx(7,j)


              l = l + 1
              enddo
          enddo
      enddo

  200 format(i4,1x,i3,3(1x,i2),3(1x,f8.3),1x,i7,69(1x,f8.2))
  300 format(A)
  100 close(7)
      close(11)
      end
