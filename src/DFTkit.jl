module dftkit
      using Shell
      using DelimitedFiles

      struct pre

      end

      struct post

      end

      struct sumr

      end

      function check_kpoint(scf_file,nk)
            scf = scf_file
            Shell.run("cp $scf tmp_scf")

            data = []
            for k in nk
                  pos = -1
                  f = open("tmp_scf") do file
                        for (iline,line) in enumerate(eachline(file))
                              if occursin("K_POINTS",line) == true
                                    pos = iline
                              end
                              if iline == pos + 1
                                    kline = " $k  $k  $k  0 0 0"
                                    cmd = "sed -i 's/$line/$kline/g' tmp_scf"
                                    Shell.run(cmd)
                                    break
                              end
                        end
                  end

                  case = " $k  $k  $k  0 0 0"
                  @info "Running for $case"
                  cmd = "mpirun -np 4 pw.x -in tmp_scf > k_analysis.out"
                  Shell.run(cmd)

                  prefix = "k_analysis.out"
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

      function get_dos(prefix,tmp)
            fildos = """&DOS
 prefix=\"$prefix\",
 outdir=\"$tmp\",
 fildos=\"tmp.dos\"/
"""

            open("fildos.in","w") do io
                  write(io,fildos)
            end
            Shell.run("dos.x -in fildos.in")

            x = readdlm("tmp.dos")
            if typeof(x[2,4]) == Float64
                  nspin = 2
            else
                  nspin = 1
            end
            @show nspin
            w = x[2:end,1]
            w = convert(Vector{Float64},w)
            if nspin == 2
                  dup = x[2:end,2]
                  dup = convert(Vector{Float64},dup)
                  ddw = x[2:end,3]
                  ddw = convert(Vector{Float64},ddw)
                  cum_occ = x[2:end,4]
                  cum_occ = convert(Vector{Float64},cum_occ)

                  Shell.run("rm -rf fildos.in")
                  Shell.run("rm -rf tmp.dos")

                  return w,dup,ddw,cum_occ
            elseif nspin == 1
                  dos = x[2:end,2]
                  dos = convert(Vector{Float64},dos)
                  cum_occ = x[2:end,3]
                  cum_occ = convert(Vector{Float64},cum_occ)

                  Shell.run("rm -rf fildos.in")
                  Shell.run("rm -rf tmp.dos")

                  return w,dos,cum_occ
            end

      end

      function get_band(prefix,pwband,tmp,spin)
            if spin == "up"
                  filband = """&BANDS
       outdir=\"$tmp\",
       prefix=\"$prefix\",
       filband=\"tmp.band\",
       spin_component = 1
      /
      """
            elseif spin == "down"
                  filband = """&BANDS
       outdir=\"$tmp\",
       prefix=\"$prefix\",
       filband=\"tmp.band\",
       spin_component = 2
      /
      """
            end

            open("filband.in","w") do io
                  write(io,filband)
            end
            Shell.run("bands.x -in filband.in")

            nb = 0
            open(pwband) do file
                  for (iline,line) in enumerate(eachline(file))
                        if occursin("nbnd",line) == true
                              nb = split(line,"=")
                              nb = split(nb[2])[1]
                              nb = parse(Int64,nb)
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

            Shell.run("rm -rf tmp.band")
            Shell.run("rm -rf tmp.band.gnu")
            Shell.run("rm -rf tmp.band.rap")

            return en,bnd
      end
end
