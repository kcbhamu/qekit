module qekit
      using Shell
      using DelimitedFiles
      using HDF5

      struct pre

      end

      struct post

      end

      function check_kpoint(scf_file,nk; offset=false)
            scf = scf_file
            Shell.run("cp $scf tmp_scf")

            if offset==false
                  off = 0
            else
                  off = 1
            end

            data = []
            for k in nk
                  pos = -1
                  f = open("tmp_scf") do file
                        for (iline,line) in enumerate(eachline(file))
                              if occursin("K_POINTS",line) == true
                                    pos = iline
                              end
                              if iline == pos + 1
                                    kline = " $k  $k  $k  $off $off $off"
                                    cmd = "sed -i 's/$line/$kline/g' tmp_scf"
                                    Shell.run(cmd)
                                    break
                              end
                        end
                  end

                  prefix = "k_analysis.out"
                  @info "Running for $k  $k  $k  $off $off $off"

                  cmd = "mpirun -np 4 pw.x -in tmp_scf > $prefix"
                  Shell.run(cmd)

                  f = open(prefix) do file
                        for (iline,line) in enumerate(eachline(file))
                              if occursin("!    total energy",line) == true
                                    en = split(line,"=")
                                    en = parse(Float64,split(en[2])[1])
                                    push!(data,en)
                              end
                        end
                  end
                  Shell.run("rm -rf $prefix")
            end
            Shell.run("rm -rf tmp_scf")
            return data
      end

      function check_ecut(scf_file,ecut)
            scf = scf_file
            Shell.run("cp $scf tmp_scf")

            data = []
            for ec in ecut
                  pos = -1
                  f = open("tmp_scf") do file
                        for (iline,line) in enumerate(eachline(file))
                              if occursin("ecutwfc",line) == true
                                    cline = "ecutwfc=$ec"
                                    cmd = "sed -i 's/$line/$cline/g' tmp_scf"
                                    Shell.run(cmd)
                              elseif occursin("ecutrho",line) == true
                                    er = ec*5
                                    cline = "ecutrho=$er"
                                    cmd = "sed -i 's/$line/$cline/g' tmp_scf"
                                    Shell.run(cmd)
                              end
                        end
                  end

                  case = "ecut = $ec"
                  @info "Running for $case"
                  cmd = "mpirun -np 4 pw.x -in tmp_scf > ecut_analysis.out"
                  Shell.run(cmd)

                  prefix = "ecut_analysis.out"
                  f = open(prefix) do file
                        for (iline,line) in enumerate(eachline(file))
                              if occursin("!    total energy",line) == true
                                    en = split(line,"=")
                                    en = split(en[2])[1]
                                    en = parse(Float64,en)
                                    push!(data,en)
                              end
                        end
                  end

                  Shell.run("rm -rf $prefix")
            end
            Shell.run("rm -rf tmp_scf")
            return data

      end

      function get_dos(; prefix, tmp, erange)
            # write fildos for dos.x input
            fildos = """&DOS
                   prefix=\"$prefix\",
                   outdir=\"$tmp\",
                   Emin=$(erange[1])
                   Emax=$(erange[2])
                   fildos=\"tmp.dos\"
                  /
                  """

            open("fildos.in","w") do io
                  write(io,fildos)
            end

            # run dos.x
            Shell.run("dos.x -in fildos.in")

            # read dos data
            x = readdlm("tmp.dos")
            if typeof(x[2,4]) == Float64
                  nspin = 2
            else
                  nspin = 1
            end

            w = x[2:end,1]
            w = convert(Vector{Float64},w)
            if nspin == 2
                  dup = convert(Vector{Float64},x[2:end,2])
                  ddw = convert(Vector{Float64},x[2:end,3])
                  cum_occ = convert(Vector{Float64},x[2:end,4])

                  Shell.run("rm -rf fildos.in")
                  Shell.run("rm -rf tmp.dos")

                  h5write("data.h5","totalDOS/E (eV)", w)
                  h5write("data.h5","totalDOS/DOS_up", dup)
                  h5write("data.h5","totalDOS/DOS_dw", ddw)
                  h5write("data.h5","totalDOS/cum_DOS", cum_occ)

                  @info "total DOS data is saved in data.h5"
            elseif nspin == 1
                  dos = convert(Vector{Float64},x[2:end,2])
                  cum_occ = convert(Vector{Float64},x[2:end,3])

                  Shell.run("rm -rf fildos.in")
                  Shell.run("rm -rf tmp.dos")

                  h5write("data.h5","totalDOS/E (eV)", w)
                  h5write("data.h5","totalDOS/DOS", dos)
                  h5write("data.h5","totalDOS/cum_DOS", cum_occ)

                  @info "total DOS data is saved in data.h5"
            end

      end

      function get_band(; prefix,pwband,tmp,spin)
            # write filband input for bands.x
            filband = """&BANDS
                   outdir=\"$tmp\",
                   prefix=\"$prefix\",
                   filband=\"tmp.band\",
                   spin_component = $spin
                  /
                  """

            open("filband.in","w") do io
                  write(io,filband)
            end
            # run bands.x
            Shell.run("bands.x -in filband.in")

            # read bands data
            nb = 0
            open(pwband) do file
                  for (iline,line) in enumerate(eachline(file))
                        if occursin("nbnd",line) == true
                              nb = split(line,"=")
                              nb = parse(Int64,split(nb[2])[1])
                        end
                  end
            end

            x = readdlm("tmp.band.gnu")
            len = Int(size(x,1)/nb)
            en = x[1:len,1]
            x = x[:,2]

            bnd = zeros(Float64,len,nb)
            for i in 1:nb
                  bnd[:,i] = x[1+len*(i-1):len+len*(i-1)]
            end

            # remove cache file
            Shell.run("rm -rf filband.in")
            Shell.run("rm -rf tmp.band")
            Shell.run("rm -rf tmp.band.gnu")
            Shell.run("rm -rf tmp.band.rap")

            h5write("data.h5","band/kpath", en)
            h5write("data.h5","band/bands", bnd)

            @info "bands data is saved in data.h5"
      end

      function get_pdos(; prefix, outdir, erange)
            function output()
                  fil = readdir()

                  idx = []
                  for (ifil,dfil) in enumerate(fil)
                        if findfirst("pdos.dat",dfil) != nothing append!(idx,ifil) end
                  end
                  return fil[idx]
            end

            # orbital_species
            s = ["s"]
            p = ["pz", "px", "py"]
            d = ["dz2", "dzx", "dzy", "dx2-y2", "dxy"]

            # write filproj input for projwfc.x
            filproj = """&PROJWFC
             prefix=\"$prefix\",
             outdir=\"$outdir\",
             Emin=$(erange[1]),
             Emax=$(erange[2]),
             filpdos=\"pdos.dat\"
            /
            """

            open("filproj.in","w") do io
                  write(io,filproj)
            end

            # run projwfc.x
            Shell.run("projwfc.x -in filproj.in")

            # output file
            fil = output()

            # save file to data.h5
            for (ifil,dfil) in enumerate(fil[1:end-1])
                  #get atom name
                  if dfil[22] == ')'
                        atom = "atom_"*dfil[19:22]
                        wfc = dfil[28]
                        orb = dfil[30]
                        data = readdlm(dfil)
                  else
                        atom = "atom_"*dfil[19:23]
                        wfc = dfil[29]
                        orb = dfil[31]
                        data = readdlm(dfil)
                  end

                  # write energy mesh to file
                  if ifil == 1
                        mesh = float.(data[2:end,1])
                        h5write("data.h5", "partialDOS/E(eV)", mesh)
                  end

                  # write pdos/ldos data
                  # number of dos data = nspin * (total_dos + orbital_species)
                  # ob = s/p/d
                  if orb == 's'
                        ob = s
                  elseif orb == 'p'
                        ob = p
                  elseif orb == 'd'
                        ob = d
                  end

                  # write for each atom
                  # size is nspin*(Ldos + norb)
                  # shift +1 because the first column is mesh
                  # start from 2 to size
                  for i in 2:( 2*(1+length(ob)) + 1 )
                        dos = float.(data[2:end,i])

                        # the first two is for Local DOS data for orbital
                        if i == 2
                              path = "partialDOS/"*atom*"/"*wfc*orb*"/LDOSup"
                              h5write("data.h5", path, dos)
                        elseif i == 3
                              path = "partialDOS/"*atom*"/"*wfc*orb*"/LDOSdw"
                              h5write("data.h5", path, dos)
                        # else is partial DOS for each orbital species
                        else
                              iob = floor(Int,i/2) - 1
                              if i % 2 == 0
                                    path = "partialDOS/"*atom*"/"*wfc*orb*"/"*ob[iob]*"_up"
                                    h5write("data.h5", path, dos)
                              else
                                    path = "partialDOS/"*atom*"/"*wfc*orb*"/"*ob[iob]*"_dw"
                                    h5write("data.h5", path, dos)
                              end
                        end
                  end #end for
            end #end for

            Shell.run("rm -rf filproj.in")
            Shell.run("rm -rf lowdin.txt")
            Shell.run("rm -rf pdos.dat.*")

            @info "Partial DOS data is saved in data.h5"
      end

end #end qekit module
