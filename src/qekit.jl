module qekit
    using Shell
    using HDF5
    using Suppressor
    using DelimitedFiles

    import Base: show

    struct qe
        out    :: Dict{String,String}
        prefix :: String
        nspin
        fermi  :: Float64
        energy :: Float64
        magnet
        charge
    end

    struct wannier90
    end

    function qe_config(prefix; scf_out)
        # Read Outdir in prefix.scf
        scf_in = prefix * ".scf"
        outdir = split(read(`grep "outdir" $scf_in`,String),"=")[2]
        outdir = split(outdir,"\"")[2]

        # Read Fermi energy
        fermi = split(read(`grep "the Fermi energy is" $scf_out`,String))[5]
        fermi = parse(Float64,fermi)

        # Read Total energy
        energy = split(read(`grep "!    total energy" $scf_out`,String),"=")
        energy = parse(Float64,split(energy[2])[1])

        # Read nspin
        try
            nspn = split(read(`grep "nspin" $scf_in`,String),"=")[2]
            if parse(Float64,nspn) == 1
                global nspn = [1]
            else
                global nspn = [1,2]
            end
        catch
            global nspn = [1]
        end

        # Read Magnet # TO DO

        # Read Charge # TO DO

        # Construct struct
        out = Dict("outdir" => outdir,
                   "scf_out" => scf_out)
        return qe(out,prefix,nspn,fermi,energy,Dict(),Dict())
    end

    #TO DO
    function show(io::IO,crystal::qe)
    end

    getMagnet(crystal::qe) = crystal.magnet
    getFermi(crystal::qe) = crystal.fermi
    getEnergy(crystal::qe) = crystal.energy
    getCharge(crysta::qe) = crystal.charge

    function pw_band(crystal::qe)
        prefix = crystal.prefix
        outdir = crystal.out["outdir"]
        nspin  = crystal.nspin

        # local variable to be used in for loop
        local kp,band
        # Start iteration for every spin up and down
        for spin in nspin
            if typeof(spin) == String
                spin = "\"none\""
            end
            filband = """&BANDS
                     outdir         = '$outdir'
                     prefix         = '$prefix'
                     filband        = 'tmp.band'
                     spin_component = $spin
                    /
                    """
            open("filband.in","w") do io
                write(io,filband)
            end

            # Run bands.x
            @info "Collecting the bands"
            @suppress run(`bands.x -in filband.in`)

            # Read bands data
            fband = readdlm("tmp.band.gnu")

            # Read number of bands
            if spin == 1
                nb = split(read(`grep "nbnd=" tmp.band`,String),"=")[2]
                nb = parse(Int64,split(nb,",")[1])
            end

            # Separate kpath and bands data
            kp  = fband[1:Int(size(fband,1)/nb),1]
            bands = fband[:,2]

            # Write bands to column vectors
            len = Int(size(fband,1)/nb)
            bnd = zeros(Float64,len,nb)
            for i in 1:nb
                bnd[:,i] = bands[1+len*(i-1):len+len*(i-1)]
            end

            # remove cache file
            run(`rm -rf filband.in tmp.band tmp.band.gnu tmp.band.rap`)

            # Concatinate the band if nspin = 2
            if spin == 1
                band = bnd
            elseif spin == 2
                dummy = zeros(Float64,2,size(band)...)
                dummy[1,:,:] = band
                dummy[2,:,:] = bnd
                band = dummy
            end
        end

        h5write("data.h5","pw_band/kpath", kp)
        h5write("data.h5","pw_band/bands", band)

        @info "bands data is saved in data.h5"
    end

    function pw_dos(crystal::qe; erange)
        prefix = crystal.prefix
        outdir = crystal.out["outdir"]
        nspin  = crystal.nspin

        # input file for dos.x
        fildos = """&DOS
                   prefix = '$prefix',
                   outdir = '$outdir',
                   Emin   = $(erange[1])
                   Emax   = $(erange[2])
                   fildos = 'tmp.dos'
                  /
                  """

        # write fildos
        open("fildos.in","w") do io
            write(io,fildos)
        end

        # run dos.x
        @info "Collecting density of states"
        @suppress run(`dos.x -in fildos.in`)

        # read dos data
        fdos = readdlm("tmp.dos")
        # get energy bins
        w = convert(Vector{Float64},fdos[2:end,1])
        # get dos data for spin up
        dos = convert(Vector{Float64},fdos[2:end,2])
        ## for spin down, and concatinate them
        tmp_dos = convert(Vector{Float64},fdos[2:end,3])
        dos = hcat(dos,tmp_dos)
        cum_occ = convert(Vector{Float64},fdos[2:end,4])

        # remove cache file
        run(`rm -rf fildos.in tmp.dos.in tmp.dos`)

        # write data to hdf5 file
        h5write("data.h5", "pw_dos/energies", w)
        h5write("data.h5", "pw_dos/dos", dos)
        h5write("data.h5", "pw_dos/cum_occ", cum_occ)
        @info "density of state of the system is saved in data.h5"
    end

    function pw_pdos(crystal::qe; erange)
        prefix = crystal.prefix
        outdir = crystal.out["outdir"]
        nspin  = crystal.nspin

        # function to detect list of all outputs file
        # that is produced by projwfc.x
        function outputs()
            fil = readdir() # list of all files
            idx = [] # empty array, for pdos index file
            for (ifil,dfil) in enumerate(fil)
                # if there are any file that contain names pdos.dat,
                # put the file into list
                if findfirst("pdos.dat",dfil) != nothing append!(idx,ifil) end
            end
            return fil[idx]
        end

        # list of all orbital species
        # TO DO for f-band
        s = ["s"]
        p = ["pz", "px", "py"]
        d = ["dz2", "dzx", "dzy", "dx2-y2", "dxy"]

        # write input file for projwfc.x
        filproj = """&PROJWFC
             prefix = '$prefix',
             outdir = '$outdir',
             Emin   = $(erange[1]),
             Emax   = $(erange[2]),
             filpdos= 'pdos.dat'
            /
            """
        open("filproj.in","w") do io
            write(io,filproj)
        end

        # run projwfc.x
        @info "Collecting projected dos"
        @suppress run(`projwfc.x -in filproj.in`)

        # list of all output files
        fil = outputs()

        # get the data for every proj dos file, except the last one
        for (ifil,dfil) in enumerate(fil[1:end-1])
            # get atom identity (name, wavefunction number, and orbital)
            ## if the atom has two words, ex : Sr, Ti, etc
            ## then there are shift on the index
            shift = 0
            if dfil[22] != ')' shift = 1 end
            atom = "atom_"*dfil[19:22+shift]
            wfc = dfil[28+shift]
            orb = dfil[30+shift]

            # convert orb to orbital species data
            if (orb == 's') vorb = s end
            if (orb == 'p') vorb = p end
            if (orb == 'd') vorb = d end

            # read the file
            data = readdlm(dfil)

            # write energy bins to file
            # only once
            if ifil == 1
                energies = float.(data[2:end,1])
                h5write("data.h5", "pdos/energies", energies)
            end

            # write pdos and ldos data to file, for each atom
            # number of total dos data for each file is (ldos + orbital_species)
            ndata_start = 2 # start from 2, because 1 is energy bins
            ndata_end   = (1 + length(vorb))
            for ildos in ndata_start:2:(2*ndata_end + 1) # jump every 2, because spins
                dosup = float.(data[2:end,ildos])
                dosdw = float.(data[2:end,ildos+1])
                dos   = hcat(dosup,dosdw)

                # the first one is for ldos
                if ildos == 2
                    path = joinpath("pdos",atom,wfc*orb,"ldos")
                else # else for pdos for each orbital
                    # after ldos, run for every orbital species
                    iorb = floor(Int,ildos/2) - 1
                    path = joinpath("pdos",atom,wfc*orb,vorb[iorb])
                end
                # read the data to file
                h5write("data.h5", path, dos)
            end
        end

        # Remove cache files
        run(`rm -rf filproj.in lowdin.txt`)
        Shell.run("rm -rf pdos.dat.*")

        @info "pdos data is saved in data.h5"
    end

    function get_w90Input(crystal::qe)
    end

    function w90_band(w90::wannier90)
    end

    function w90_pdos(w90::wannier90)
    end

    function w90_boltz(w90::wannier90)
    end

    function w90_berry(w90::wannier90)
    end

    function qekit_info()
    end

    function w90_lib()
    end

    function pw_lib()
    end

end
