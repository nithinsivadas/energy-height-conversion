      program testRead

      ! The program takes two arguments, specifying the input
      ! filename and the output filename (including paths.)
      ! The output file is written in the current working directory
      ! unless a path is included as part of the filename.

      integer*2 j
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

      do rec = 1,5400
          call readArcSub(7,rec,eof)
          if(eof.eq.1) goto 100
          do j=1,4
              hr = ihd(4,j)/3600000
              min = (ihd(4,j) - hr*3600000)/60000
              sec = (ihd(4,j) - hr*3600000 - min*60000)/1000
              write(11,200),ihd(2,j),ihd(3,j),hr,min,sec,
     +                  head(2,j),head(3,j),ihd(6,j)
          enddo
          write(11,300)"mep0:"
          do j=1,16
              write(11,220)j,mep0(1,j),mep0(2,j),mep0(3,j),mep0(4,j),
     +          mep0(5,j),mep0(6,j),mep0(7,j),mep0(8,j),mep0(9,j)
          enddo
          write(11,300)"mep90:"
          do j=1,16
              write(11,220)j,mep90(1,j),mep90(2,j),mep90(3,j),
     +          mep90(4,j),mep90(5,j),mep90(6,j),mep90(7,j),mep90(8,j),
     +          mep90(9,j)
          enddo
          write(11,300)"mepOmni:"
          do j=1,16
             write(11,220)j,mepOmni(1,j),mepOmni(2,j),mepOmni(3,j),
     +       mepOmni(4,j)
          enddo
          write(11,300)"ted0:"
          do j=1,16
              write(11,220)j,ted0(1,j),ted0(2,j),ted0(3,j),ted0(4,j),
     +          ted0(5,j),ted0(6,j),ted0(7,j),ted0(8,j)
          enddo
          write(11,300)"ted30:"
          do j=1,16
              write(11,220)j,ted30(1,j),ted30(2,j),ted30(3,j),
     +          ted30(4,j),ted30(5,j),ted30(6,j),ted30(7,j),
     +          ted30(8,j)
          enddo
          write(11,300)"ted0s:"
          do j=1,4
              write(11,220)j,ted0s(1,j),ted0s(2,j),ted0s(3,j),
     +          ted0s(4,j),ted0s(5,j),ted0s(6,j),ted0s(7,j),ted0s(8,j)
          enddo
          write(11,300)"ted30s:"
          do j=1,4
              write(11,220)j,ted30s(1,j),ted30s(2,j),ted30s(3,j),
     +          ted30s(4,j),ted30s(5,j),ted30s(6,j),ted30s(7,j),
     +          ted30s(8,j)
          enddo
          write(11,300)"tedback:"
          do j=1,4
              write(11,220)j,tedback(1,j),tedback(2,j)
          enddo
          write(11,300)"tedfx:"
          do j=1,16
              write(11,220)j,tedfx(1,j),tedfx(2,j),tedfx(3,j),
     +           tedfx(4,j),tedfx(5,j),tedfx(6,j),tedfx(7,j)
          enddo

          write(11,300)"----------------------------------------------"

      enddo

  200 format(i4,1x,i3,3(1x,i2),1x,f8.3,1x,f8.3,1x,i7)
  220 format(i4,9(1x,f8.2))
  300 format(A)
  100 close(7)
      close(11)
      end
