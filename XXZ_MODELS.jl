
using LinearAlgebra
using Statistics
using StatsBase
BLAS.set_num_threads(4)
using ITensors, ITensorMPS
include("Models.jl")

using .Models
let
  # Define number of spins and create spin-one indices
  N = parse(Int,ARGS[1])
  W = parse(Int,ARGS[2])
  h = parse(Float64,ARGS[3])
  hf = parse(Float64,ARGS[4])
  dh = parse(Int,ARGS[5])
  NS = parse(Int,ARGS[6])
  NR = parse(Int,ARGS[7])
  ang = parse(Int,ARGS[8])
  model = ARGS[9]
  NC = 0.2
  sites = siteinds("S=1/2", N*W)
  # Initialize a random MPS
  psi0 = random_mps(sites)
  Nh=Int((hf-h)/dh)
  Nang=10.0
  
  var=.0
  folder="RESULTS"*"/"*"RUN_QIM_DMRG_L" * ARGS[1] * "W" * ARGS[2]* "NS" * ARGS[5] * "G" * ARGS[3] * "GF" * ARGS[4] * "NR" * ARGS[6] * "ANG" * "$ang" * model

  mkpath(folder)

  filename_2 = folder*"/"*"DATAM5L" * ARGS[1] * "W" * ARGS[2]* "NS" * ARGS[6] * model * "MPS" * ".txt" # Output file optimization

  
  open(filename_2,"a") do file
              write(file,"G","\t","E","\t","E_var","\n")
  end

        	      
   
  for iter in 0:Nh

      h_model=h+iter*dh

      En=zeros(NR)
      
      println("L=$N","W=$W","G=$h_model","NS=$NS ",model)
      os=rotated_over_y_quantum_heisenberg(N,h_model,float(ang*pi)/Nang)
      H = MPO(os, sites)
      
      
      for sample_iter in 1:NR
      
      	  filename_opt = folder*"/"*"OPTDATAM5L" * ARGS[1] * "W" * ARGS[2] * "NS" * ARGS[6] * "MPSG" * string(Int(h_model)) * ".txt" * string(sample_iter)  # Output file to store configurations
	  filename = folder*"/"*"DATAM5L" * ARGS[1] * "W" * ARGS[2] * "NS" * ARGS[6] * "MPSG" * string(Int(h_model)) * ".txt" * string(sample_iter)  # Output file to store configurations
	  
	  # Run DMRG to find the ground state

	  nsweeps = 10
      	  maxdim = [64,64,64,128,256,256,256,400,400,512,1024,1024,1024,1024,1024]
      	  cutoff = 1E-10
	  energy, psi = DRIVER(filename_opt,H,psi0,N,nsweeps,cutoff)

	  # Update the wavefunction

	  psi0=psi
	  
	  # Compute the variance
	  
	  H2 = inner(H,psi,H,psi)
	  
	  var = H2-energy^2
	  
	  
	  En[sample_iter]=energy
      	  println("Final energy = $energy","Final var =$var")
      	  

      	  # Call the function to sample and save to the file
      	  sample_mps_to_file(psi, filename, NS)
      	  	nothing
		
	  NE=Int( floor(NC*length(psi)))
	  indices_to_modify=rand(1:length(psi0),NE)
	  for K in indices_to_modify
    	      tensor_data = ITensors.data(psi0[K])  # Access raw data
    	      tensor_data .+= randn(eltype(tensor_data), size(tensor_data))
    	      psi0[K] = ITensor(tensor_data, inds(psi0[K])...)  # Reconstruct the ITensor with modified data
	  end

	  
      end
      
      E_mean=mean(En)
      E_var=real(var)
      
      
      open(filename_2,"a") do file

              write(file,"$h_model","\t","$E_mean","\t","$E_var","\n")
      	  end
  end	
  
end
