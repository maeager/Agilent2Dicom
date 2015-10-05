#!/usr/bin/env python
"""ReconFID is used to read and process Agilent FID image files

  This is an attempt to re-create Xrecon's EPI and FSE resconstruction

   - Michael Eager (michael.eager@monash.edu)
"""
"""
  Copyright (C) 2014 Michael Eager

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

#-----------------------------*/
#---- Some sequence modes ----*/
#-----------------------------*/
# For 1D */
IM1D =   100
# For 2D multislice 200 < seqmode < 300 */
IM2D =      200
IM2DCC =    201  # seqcon = "*ccnn" */
IM2DCS =    202  # seqcon = "*csnn" */
IM2DSC =    203  # seqcon = "*scnn" */
IM2DSS =    204  # seqcon = "*ssnn" */
IM2DCCFSE = 205  # apptype = "im2Dfse", seqcon = "nccnn" 
IM2DCSFSE = 206  # apptype = "im2Dfse", seqcon = "ncsnn" 
IM2DSCFSE = 207  # apptype = "im2Dfse", seqcon = "nscnn" 
IM2DSSFSE = 208  #* apptype = "im2Dfse", seqcon = "nssnn" 
IM2DEPI =   209  # apptype = "im2Depi" 
# For 3D 300 < seqmode < 400 
IM3D =      300
IM3DCC =    301  # seqcon = "**ccn" 
IM3DCS =    302  # seqcon = "**csn" 
IM3DSC =    303  # seqcon = "**scn" 
IM3DSS =    304  # seqcon = "**ssn" 
IM3DCFSE =  305  # apptype = "im3Dfse", seqcon = "ncccn" 
IM3DSFSE =  306  # apptype = "im3Dfse", seqcon = "nccsn" 
# For 4D 
IM4D   =   400

# Knoisefraction: Fraction of k-space FOV to use to sample noise #
Knoisefraction=0.05

# IMnoisefraction: Fraction of image space FOV to use to sample noise #
IMnoisefraction=0.05
procpar=dict()
fid_header=dict()


OFF = 0
POINTWISE = 1
TRIPLE = 2
SCALED_TRIPLE = 3
  

def sliceorder(dim,par):

  tic()
  #sprocpar['ns']
  # Set dimorder=-1 if par is not found #
  standardorder=True
  if par not in procpar.keys():
      dimorder[0]=-1
  else:
    # Take a copy of the starting slice order #
    pss = procpar[par]

    # Sort slice positions into ascending order #
    pss_sorted = np.sort(pss)
    
    # Fill dimorder with an index of where successive slices are #

    for i in xrange(0,dim):
        for j in xrange(0,dim):
            if pss[i] == pss[j]:
                dimorder[i]=j
                break
    # If any slices have equal positions their dimorder index will be equal #
    # Adjust any duplicate values appropriately #
   for i in xrange(0,dim):
       k=1
       for j in xrange(i+1,dim):
           if dimorder[j]==dimorder[i]:
                dimorder[j]+=k
                k+=1

    # Set dimorder=-1 if we simply have standard order #
    for i in xrange(0,dim):
        if dimorder[i] != i:
            standardorder=False
            

    if standardorder:
        del dimorder
        dimorder = -1

    if args.verbose:
        print toc()
        print pss,dimorder
  return dimorder



def phaseorder(procpar, views, dim, par):
    
    tic()

    # Set dimorder=-1 if par is not found or has wrong size #
    standardorder=True
    if par in procpar.keys():
        porder = procpar[pars]
        # Initialize the values to be out of range #
        dimorder = np.ones(len(porder))*dim+1
       
    # Set dimorder=-1 if par has too few values #
    if dim < views:
        dimorder[0]=-1
    else:
        # Fill dimorder with an index of where successive phase encodes are #

      # If there's more values than views set the standard minimum multiplier #
      if (dim>views):
          mindim = -views/2+1
      else:
        # Find the minimum value (phase encode multipliers can run from -nv/2 to +nv/2-1 or -nv/2+1 to +nv/2) #
        mindim=0
        for i in xrange(0,dim):
            if (porder[i]<mindim):
                mindim=porder[i]

      for i in xrange(0,dim):
          for j in xrange(0,dim):
              if porder[i] == j+mindim:
                  dimorder[i]=j
                  break

      # Set dimorder=-1 if we simply have standard order #
      for i in xrange(0,dim):
          if dimorder[i] != i:
              standardorder=False
              break


      if standardorder:
          clear dimorder
          dimorder[0]=-1


  if args.verbose:
    print " Phase ordering took %f secs " % (toc())
    print dimorder

  return dimorder


def setnvolsEPI():


  # Set nuber of "volumes" #
  fid_header["volumes"]=fid_header["nblocks"]/fid_header["nr"]


  if (procpar["allvolumes"] == "n"):  # Don't process all volumes #
    fid_header["startvol"]=procpar["startvol"]
    fid_header["endvol"]=procpar["endvol"]
    if fid_header["startvol"] > fid_header["endvol"]:      # Swap them #
      fid_header["vol"]=fid_header["endvol"]; fid_header["endvol"]=fid_header["startvol"];
      fid_header["startvol"]=fid_header["vol"]; fid_header["vol"]=0

    if (fid_header["startvol"] < 0):
      fid_header["startvol"]=0
    if (fid_header["endvol"] > fid_header["volumes"]):
      fid_header["endvol"]=fid_header["volumes"]
  else: # Process all volumes #
      fid_header["startvol"]=0
      fid_header["endvol"]=fid_header["volumes"]


def defaultEPI(procpar,fid_header):
  DISCARD=-1
  refEPI=False
  refSGE=False
  getscaleref=False
  #double oversample
  #struct data ref1,ref2,ref3,ref4,ref5,ref6
  #struct segscale scale1,scale2




  fid_header["volumes"]=fid_header["nblocks"]/fid_header["nr"]
  fid_header["startvol"]=0
  fid_header["endvol"]=fid_header["volumes"]


  fid_header["nv"] = procpar["nphase"];    # Set fid_header["nv"] for dimorder2D #
  fid_header["pssorder"] = sliceorder(d,fid_header["ns"],"pss"); # Fill pssorder with the slice order #
  fid_header["dim2order"] = phaseorder(procpar,fid_header["nv"],fid_header["nv"],"par_does_not_exist"); # Set dim2order=-1 for sequential phase encode order #
  fid_header["dim3order"]=phaseorder(d,fid_header["nv"],fid_header["nv"],"sgepelist"); # Fill dim3order with the standard gradient echo phase encode order #
  fid_header['nv2']=fid_header["nv"];                       # Set procpar['nv2'] for dim3order #
  fid_header['nv']=procpar["nseg"];      # Use fid_header["nv"] for the number of shots #

  # Set EPI correction scheme #
  if (procpar["epi_pc"] == "POINTWISE"):
    epi_pc=POINTWISE
  elif (procpar["epi_pc"] == "TRIPLE"):
    epi_pc=TRIPLE
  elif (procpar["epi_pc"] == "SCALED_TRIPLE"):
    epi_pc=SCALED_TRIPLE
  else:
    epi_pc=OFF

  # Check whether to output or discard reference scans #
  if (procpar["imRF"] == "y"):
    refEPI=True
  if (procpar["imSGE"] == "y"):
    refSGE=True

  # Set reference data #
  if (epi_pc > OFF) :        # Pointwise or triple reference phase correction #
    ref1 = fid_header.copy()
    ref3 = fid_header.copy()   # ref3 used in pointwise phase correction of inverted
                      #          reference scans if reference output is selected 

  if (epi_pc > POINTWISE) :  # Triple reference phase corrections #
    ref2= fid_header.copy()
    ref4= fid_header.copy()

  if (epi_pc > TRIPLE) : # Scaled triple reference phase correction #
     ref5 = fid_header.copy()
     ref6 = fid_header.copy()


  # Set default of no compressed segment scaling just in case it's never set #
  scale1.data=False
  scale2.data=False

  # Loop over data blocks #
  for fid_header["block"] in xrange(0,fid_header["nblocks"]):
    fid_header["outvol"]=0 # Initialise output data volume (that will not include reference scans) #

    for fid_header["vol"] in xrange (0,fid_header["volumes"]): # Loop over "volumes" #

      setoutvolEPI(d)               # Set output data volume #

      if (fid_header["outvol"]>fid_header["endvol"]):
        break # Break if last requested volume has been processed #

      getblockEPI(d,fid_header["vol"],NDCC);     # Get block without applying dbh.lvl and dbh.tlt #
      zeromax(d);                     # Zero max structure & coordinates of maximum #
      zeronoise(d);                   # Zero values in noise structure #

      if epi_pc ==  OFF :
        print "  Correction: OFF \n"
      elif epi_pc == POINTWISE:
        print "  Correction: POINTWISE \n"
      elif epi_pc == TRIPLE:
        print "  Correction: TRIPLE \n"
      elif epi_pc == SCALED_TRIPLE:
        print "  Correction: SCALED_TRIPLE \n"

      # Process data as directed by image parameter #
      if procpar["image"][fid_header["curvol"]] ==  0:   # Reference, phase-encode #
        setblockEPI(d);                     # Set block for 2D data (fid_header["nv"] > 1) and navigators #
        #if (refEPI) w2Dfdfs(d,VJ,FLT32,fid_header["vol"]); # Output raw data for the volume, requested #
        if (epi_pc > OFF) :                 # If there is phase correction #
          ftnpEPI(d)                       # FT along readout dimension #
          clear2Ddata(ref1)               # Clear ref1 #
          copy2Ddata(d,ref1)              # Copy to ref1 #
          ref1.datamode=EPIREF             # Flag as EPIREF data #
        else:                            # else there is no phase correction #
          ref1.datamode=None               # Flag as no data #
          if (refEPI):
            ftnpEPI(d)           # FT along readout dimension #

        scale1=setsegscale(d)             # Set scaling for compressed segments #
        segscale(d,&scale1);                # Scale compressed segments #
        if (refEPI):                         # If reference output is requested #
          w2Dfdfs(d,VJ,FLT32,fid_header["vol"]);       # Output data for the volume #
        else:
          w2Dfdfs(d,VJ,FLT32,DISCARD);   # Else use DISCARD to flag skip the volume #
        wnifti(d,VJ,FLT32,DISCARD)
          
      elif procpar["image"][fid_header["curvol"]] == -1:  # Inverted Readout Reference, phase-encode #
#ifdef DEBUG
        print "  Processing reference -1 data ...\n"
        
#endif
        setblockEPI(d);                     # Set block for 2D data (fid_header["nv"] > 1) and navigators #
        if (refEPI):
          w2Dfdfs(d,VJ,FLT32,fid_header["vol"]); # Output raw data for the volume, requested #
        ftnpEPI(d);                         # FT along readout dimension #
        segscale(d,&scale2);                # Scale compressed segments #
        if (epi_pc > POINTWISE) :           # if triple or scaled triple reference phase correction #
          clear2Ddata(&ref2);               # Clear ref2 #
          copy2Ddata(d,&ref2);              # Copy to ref2 #
          ref2.datamode=EPIREF;             # Flag ref2 as EPIREF data #
          if (ref3.datamode == EPIREF) :    # if there is ref3 reference data #
            phaseEPIref(&ref2,&ref3,&ref4); # Phase correct ref2 data using ref3 and put result in ref4 #
            #              analyseEPInav(&ref4);           # Analyse the navigators #
            stripEPInav(&ref4);             # Strip the navigator scans #
            ftnvEPI(&ref4);                 # FT along phase encode dimension #
            revreadEPI(&ref4);              # Reverse the data in readout dimension #
            getscaleref=True;               # Flag to store the next regular image for scaling in SCALED_TRIPLE #


          if (refEPI) :                       # if reference output is requested #
            if (ref3.datamode == EPIREF) :    # if there is ref3 reference data #
              phaseEPI(d,&ref3);              # Phase correct with the reference #

            navcorrEPI(d);                    # Phase correct with the navigator #
            stripEPInav(d);                   # Strip the navigator scans #
            ftnvEPI(d);                       # FT along phase encode dimension #
            revreadEPI(d);                    # Reverse the data in readout dimension #
            w2Dfdfs(d,VJ,FLT32,fid_header["vol"]);       # Output data for the volume #
          else:
            w2Dfdfs(d,VJ,FLT32,DISCARD);   # Use DISCARD to flag skip the volume #
          wnifti(d,VJ,FLT32,DISCARD)
          
      elif    procpar["image"][fid_header["curvol"]] ==  -2:  # Inverted Readout Reference, phase-encode #
#ifdef DEBUG
        fprintf(stdout,"  Processing reference -2 data ...\n")
        fflush(stdout)
#endif
        setblockEPI(d);                     # Set block for 2D data (fid_header["nv"] > 1) and navigators #
        if refEPI:
          w2Dfdfs(d,VJ,FLT32,fid_header["vol"])
        ftnpEPI(d);                         # FT along readout dimension #
        setsegscale(d,&scale2);             # Set scaling for compressed segments #
        segscale(d,&scale2);                # Scale compressed segments #
        if (epi_pc > POINTWISE):           # if old triple or triple reference phase correction #
          clear2Ddata(&ref3);               # Clear ref3 #
          copy2Ddata(d,&ref3);              # Copy to ref3 #
          ref3.datamode=EPIREF;             # Flag ref3 as EPIREF data #
          if (ref2.datamode == EPIREF) :    # if there is ref2 reference data #
            phaseEPIref(&ref2,&ref3,&ref4); # Phase correct ref2 data using ref3 and put result in ref4 #
            #              analyseEPInav(&ref4);           # Analyse the navigators #
            stripEPInav(&ref4);             # Strip the navigator scans #
            ftnvEPI(&ref4);                 # FT along phase encode dimension #
            revreadEPI(&ref4);              # Reverse the data in readout dimension #
            getscaleref=True;               # Flag to store the next regular image for scaling in SCALED_TRIPLE #


          if (refEPI) :                       # if reference output is requested #
            if (epi_pc == POINTWISE) :        # if pointwise reference phase correction #
              clear2Ddata(&ref3);             # Clear ref3 #
              copy2Ddata(d,&ref3);            # Copy to ref3 #
              ref3.datamode=EPIREF;           # Flag ref3 as EPIREF data #

            revreadEPI(d);                    # Reverse the data in readout dimension #
            w2Dfdfs(d,VJ,FLT32,fid_header["vol"]);       # Output data for the volume #

          else:
            w2Dfdfs(d,VJ,FLT32,DISCARD);   # Use DISCARD to flag skip the volume #
          wnifti(d,VJ,FLT32,DISCARD)
          break
      elif procpar["image"][fid_header["curvol"]] ==  1:   # Regular image #
#ifdef DEBUG
        fprintf(stdout,"  Processing image 1 data ...\n")
        fflush(stdout)
#endif
        setblockEPI(d);                     # Set block for 2D data (fid_header["nv"] > 1) and navigators #
        if fid_header["outvol"]>=fid_header["startvol"]:
          w2Dfdfs(d,VJ,FLT32,fid_header["vol"]);       # Output raw data for the volume, requested #
        if epi_pc== SCALED_TRIPLE:               # Scaled triple reference phase correction #
          if (getscaleref):             # If the scale reference has just been acquired #
            clear2Ddata(&ref5);           # Clear ref5 #
            copy2Ddata(d,&ref5);          # Copy to ref5 #
            ref5.datamode=EPIREF;         # Flag ref5 as EPIREF data #
            prepEPIref(&ref5,&ref1);      # Prepare ref5 data so it can be used to scale ref4 data #

              
            
              

          ftnpEPI(d);                         # FT along readout dimension #
          segscale(d,&scale1);                # Scale compressed segments #
          phaseEPI(d,&ref1);                  # Phase correct with the reference #
          navcorrEPI(d);                      # Phase correct with the navigator #
#          analyseEPInav(d);                   # Analyse the navigators #
          stripEPInav(d);                     # Strip the navigator scans #
          ftnvEPI(d);                         # FT along phase encode dimension #
          if epi_pc == TRIPLE:                      # Triple reference phase correction #
            addEPIref(d,&ref4);             # Add ref4 data to cancel N/2 ghost #
              
          elif epi_pc == SCALED_TRIPLE:               # Scaled triple reference phase correction #
            if (getscaleref) :              # If the scale reference has just been acquired #
              addEPIref(d,&ref4);           # Add ref4 data to cancel N/2 ghost #
              getscaleref=False;            # Flag that scale reference has been acquired #
            else:
              addscaledEPIref(d,&ref4,&ref5); # Scale ref4 data according to d/ref5, add to d #


          if fid_header["outvol"]>=fid_header["startvol"]:
            phasedata2D(d,VJ);                # Phase data if required #
            w2Dfdfs(d,VJ,FLT32,fid_header["vol"]);       # Write 2D fdf data from volume #
          wnifti(d,VJ,FLT32,fid_header["vol"])

          

      else: # Reference Standard Gradient Echo #
          setblockSGE(d)
          if (refSGE):
            w2Dfdfs(d,VJ,FLT32,fid_header["vol"]); # Output raw data for the volume, requested #
          shiftdata2D(d,STD);                 # Shift FID data for fft #
          zeronoise(d);                       # Zero any noise measurement #
          equalizenoise(d,STD);               # Scale for equal noise in all receivers #
          phaseramp2D(d,READ);                # Phase ramp the data to correct for readout offset pro #
          phaseramp2D(d,PHASE);               # Phase ramp the data to correct for phase encode offset ppe #
          weightdata2D(d,STD);                # Weight data using standard VnmrJ parameters #
          zerofill2D(d,STD);                  # Zero fill data using standard VnmrJ parameters #
          fft2D(d,STD);                       # 2D fft #
          phasedata2D(d,VJ);                  # Phase data if required #
          shiftdata2D(d,STD);                 # Shift data to get images #
          oversample=procpar["oversample"]; # Check to see if there is oversampling #
          if (oversample>1):
            zoomEPI(d):       # If oversampled, to get the requested FOV #
          if (refSGE):                         # If standard gradient echo reference output is requested #
            w2Dfdfs(d,VJ,FLT32,fid_header["vol"]);       # Output data for the volume #
          else w2Dfdfs(d,VJ,FLT32,DISCARD);   # Else use DISCARD to flag skip the volume #
          wnifti(d,VJ,FLT32,DISCARD)
     # end image parameter switch #

      clear2Ddata(d);               # Clear data volume from memory #
      setdim(d);                    # Reset data dimensions in case data has been zerofilled (sets fid_header["nv"]=1) #
      fid_header["nv"]=procpar["nseg"];     # Use fid_header["nv"] for the number of shots #
      procpar['nv2']=procpar["nphase"]; # Use procpar['nv2'] for number of standard gradient echo phase encodes #
      fid_header["dimstatus"][0] = None;       # Make sure ZEROFILL status is not set, setdim will be called in getblock() #
      fid_header["dimstatus"][1] = None;       # Make sure ZEROFILL status is not set, setdim will be called in getblock() #

  
  # Clear all reference data #
  if (epi_pc > OFF) :        # Pointwise or triple reference phase correction #
    clear2Dall(&ref1)
    clear2Dall(&ref3)

  if (epi_pc > POINTWISE) :  # Triple and scaled triple reference phase correction #
    clear2Dall(&ref2)
    clear2Dall(&ref4)

  if (epi_pc > TRIPLE) :     # Scaled triple reference phase correction #
    clear2Dall(&ref5)
    clear2Dall(&ref6)


  clear2Dall(d);             # Clear everything from memory #

  return  procpar, d, image_data


  
def ParseEPIP(procpar,hdr,kimage_data):
    
    # reconEPI
    if procpar['recon'] == "prescanEPI":
      procpar, d, image_data = prescanEPI(procpar,hdr,kimage_data)
    else:
      procpar, d, image_data = defaultEPI(procpar,hdr,kimage_data)



