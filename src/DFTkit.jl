module dftkit
      using Shell

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

end
