
module Models

export sample_mps_to_file, cluster_ising, quantum_ising, quantum_XYZ,process_model,O_Z,rotated_quantum_ising,quantum_heisenberg_Z,rotated_over_y_quantum_heisenberg,DRIVER

using ITensors, ITensorMPS

function sample_mps_to_file(psi::MPS, filename::String, N::Int)
    # Open a file for writing
    open(filename, "w") do file
        for _ in 1:N
            # Sample a configuration from the MPS
            config = sample!(psi)
	    # config.-= 3/2
	    # config.*=2
	    config.=round.(Int,(config.-3/2)*2)
            # Convert the configuration into a string
            # The configuration `config` is a vector of spins at each site
            config_str = join(string.(config), " ")  # Convert each spin to a string
            println(file, config_str)  # Write to file
        end
    end
    println("Sampling complete. Configurations saved to $filename")
end

#MODELS:


function cluster_ising_z(N::Int,W::Int,h::Float64)
  os = OpSum()
  for j = 0:N-3
     os += -8.0,"Sx",j+1,"Sz",(j+1)%N+1,"Sx",(j+2)%N+1
     os += h*0.004, "Sy",j+1,"Sy",(j+1)%N+1
  end
  os+=h*0.004,"Sy",N-1,"Sy",N
  return os
end


function cluster_ising_x(N::Int,W::Int,h::Float64)
  os = OpSum()
  for j = 0:N-3
     os += -8.0,"Sz",j+1,"Sx",(j+1)%N+1,"Sz",(j+2)%N+1
     os += h*0.004, "Sy",j+1,"Sy",(j+1)%N+1
  end
  os+=h*0.004,"Sy",N-1,"Sy",N
  return os
end



function cluster_ising_y(N::Int,W::Int,h::Float64)
  os = OpSum()
  for j = 0:N-3
     os += 8.0,"Sx",j+1,"Sy",(j+1)%N+1,"Sx",(j+2)%N+1
     os += h*0.004, "Sz",j+1,"Sz",(j+1)%N+1
  end
  os+=h*0.004,"Sz",N-1,"Sz",N
  return os
end




function rotated_quantum_ising(N::Int,W::Int,h::Float64,theta::Float64)
  os = OpSum()
  if W>1
   for i = 0:W-2
     for j = 0:N-1
      os += -4.0*cos(theta),"Sz",j+1+i*N, "Sz",(j+1)%N+1+i*N
      os += -4.0*cos(theta),"Sz",j+1+i*N, "Sz",j+1+((i+1)%W)*N
      os +=  4.0*sin(theta),"Sx",j+1+i*N, "Sx",(j+1)%N+1+i*N
      os +=  4.0*sin(theta),"Sx",j+1+i*N, "Sx",j+1+((i+1)%W)*N

      os += -h*0.002*cos(theta), "Sx", j+1+i*W
      os += -h*0.002*sin(theta), "Sz", j+1+i*W
      
     end
   end
  else
   for j = 0:N-1
   
     os += -4.0*cos(theta)*cos(theta),"Sz",j+1, "Sz",(j+1)%N+1

     os += +4.0*sin(theta)*cos(theta),"Sz",j+1, "Sx",(j+1)%N+1
     os += +4.0*sin(theta)*cos(theta),"Sx",j+1, "Sz",(j+1)%N+1

     os += -4.0*sin(theta)*sin(theta),"Sx",j+1, "Sx",(j+1)%N+1

     os += -h*0.002*cos(theta), "Sx",j+1
     os += -h*0.002*sin(theta), "Sz",j+1
     
   end
  end
  return os
end






function rotated_over_y_quantum_heisenberg(N::Int,h::Float64,theta::Float64)
  os = OpSum()

  delta=0.0001*h
  
  for j = 0:N-1
   
     os += 4.0*((sin(theta))^2+((cos(theta))^2)*delta),"Sz",j+1, "Sz",(j+1)%N+1
     os += 4.0*((cos(theta))^2+delta*(sin(theta))^2),"Sx",j+1, "Sx",(j+1)%N+1
     os += 4.0,"Sy",j+1, "Sy",(j+1)%N+1

     os += (1-delta)*sin(theta)*cos(theta)*4.0,"Sx",j+1, "Sz",(j+1)%N+1
     os += (1-delta)*sin(theta)*cos(theta)*4.0,"Sz",j+1, "Sx",(j+1)%N+1

  end
   
  return os
end




function quantum_heisenberg_X(N::Int,W::Int,h::Float64)
  os = OpSum()
  if W>1
   for i = 0:W-1
     for j = 0:N-1
      os += -0.04*h,"Sx",j+1+i*W, "Sx",(j+1)%N+1+i*N
      os += -0.04*h,"Sx",j+1+i*W, "Sx",j+1+((i+1)%W)*N
      os += -4.0,"Sz",j+1+i*W, "Sz",(j+1)%N+1+i*N
      os += -4.0,"Sz",j+1+i*W, "Sz",j+1+((i+1)%W)*N
      os += -4.0,"Sy",j+1+i*W, "Sy",(j+1)%N+1+i*N
      os += -4.0,"Sy",j+1+i*W, "Sy",j+1+((i+1)%W)*N
     end
   end
   
  else
   for j = 0:N-1
     os += -0.04*h,"Sx",j+1, "Sx",(j+1)%N+1
     os += -4.0,"Sz",j+1, "Sz",(j+1)%N+1
     os += -4.0,"Sy",j+1, "Sy",(j+1)%N+1
   end
  end
  return os
end

function quantum_heisenberg_Y(N::Int,W::Int,h::Float64)
  os = OpSum()
  if W>1
   for i = 0:W-1
     for j = 0:N-1
      os += -0.04*h,"Sx",j+1+i*W, "Sx",(j+1)%N+1+i*N
      os += -0.04*h,"Sx",j+1+i*W, "Sx",j+1+((i+1)%W)*N
      os += -4.0,"Sy",j+1+i*W, "Sy",(j+1)%N+1+i*N
      os += -4.0,"Sy",j+1+i*W, "Sy",j+1+((i+1)%W)*N
      os += -4.0,"Sz",j+1+i*W, "Sz",(j+1)%N+1+i*N
      os += -4.0,"Sz",j+1+i*W, "Sz",j+1+((i+1)%W)*N
     end
   end
   
  else
   for j = 0:N-1
     os += -0.04*h,"Sx",j+1, "Sx",(j+1)%N+1
     os += -4.0,"Sy",j+1, "Sy",(j+1)%N+1
     os += -4.0,"Sz",j+1, "Sz",(j+1)%N+1
   end
  end
  return os
end


function O_Z(psi::MPS, N::Int)
    psi_R=deepcopy(psi)
    aux=N-2
    aux_2=N-1
    s = siteind("S=1/2")
    # Expectation values for sigma_x, sigma_y, sigma_z on specific sites
    for i in 3:N-2
    	psi_R = orthogonalize(psi,i)
        NEW_INDEX = op("Sz",s) * psi_R[i]
    	NEW_INDEX = noprime(NEW_INDEX)
	psi_R[i] = NEW_INDEX
    end    

    for i in [1,N]
    	psi_R = orthogonalize(psi,i)
        NEW_INDEX = op("Sx",s) * psi_R[i]
    	NEW_INDEX = noprime(NEW_INDEX)
	psi_R[i] = NEW_INDEX
    end

    for i in [2,N-1]
    	psi_R = orthogonalize(psi,i)
        NEW_INDEX = op("Sy",s) * psi_R[i]
    	NEW_INDEX = noprime(NEW_INDEX)
	psi_R[i] = NEW_INDEX
    end
    
    println(siteinds(psi))
    println(siteinds(psi_R))
    result=inner(psi,psi_R)
   
    return result
    
end

function CHAIN_O_Z(psi::MPS, N::Int)
    psi_R=deepcopy(psi)
   
    s = siteind("S=1/2")
    # Expectation values for sigma_x, sigma_y, sigma_z on specific sites
    for i in 3:N-2
    	psi_R = orthogonalize(psi,i)
        NEW_INDEX = op("Sz",s) * psi_R[i]
    	NEW_INDEX = noprime(NEW_INDEX)
	psi_R[i] = NEW_INDEX
    end
    println(siteinds(psi))
    println(siteinds(psi_R))
    result=inner(psi,psi_R)
    return result
       
        
end


function DRIVER(filename::String, H::MPO, psi0::MPS, N::Int, nsweeps, cutoff)
    psi = psi0
    energy=.0
    var=.0
    maxdim = [64,64,64,128,256,256,256,400,400,512,1024,1024,1024,1024,1024]
    for iter in 1:nsweeps

        bond_dim_idx = min(iter, length(maxdim))
        
        # Call dmrg with correct keyword arguments
        energy, psi = dmrg(H, psi; nsweeps=1, maxdim=maxdim[bond_dim_idx], cutoff=cutoff)

        # Compute variance
        H2 = inner(H,psi, H,psi)
        var = H2 - energy^2

        # Append to file
        open(filename, "a") do file
            write(file, string( real(energy), "\t", real(var), "\t", maxdim[bond_dim_idx], "\n"))
        end
    end
    
    return energy, psi
end



end

	
