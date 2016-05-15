require 'bio'

def coord(hetatm) # coordinate | return vector [x y z] 
	return Bio::PDB::Utils.to_xyz(hetatm)
end

def dihedral_angle(hash, a1, a2, a3, a4)
	if hash.has_key?(a1) and hash.has_key?(a2) and hash.has_key?(a3) and hash.has_key?(a4) # atoms are available
		return Bio::PDB::Utils.dihedral_angle(
			coord(hash[a1]),
			coord(hash[a2]),
			coord(hash[a3]),
			coord(hash[a4]))
	else
		return
	end
end

def calculate_sam(file)
	pdb = Bio::PDB.new(File.open(file).read)
	pdb_id = pdb.entry_id()

	# filter "HETATM - SAM"
	sam = pdb.hetatms.delete_if { |element| element.resName != 'SAM'}

	data = Hash.new
	sam.each do |element| # create hash based on 'name' of the atom
		if not data.has_key?(element.chainID)
			data.merge!({element.chainID => {} })
		end
		data[element.chainID].merge!({element.name => element}) # merge
	end
# alpha    CB -CG -SD -C5'
# beta     CG -SD -C5'-C4'
# gamma    SD -C5'-C4'-O4'
# delta    C5'-C4'-C3'-O3'
# epsilon  
# zeta
# chi      O4'-C1'-N9 -C2'
	data.each do |key, value|
		alpha = dihedral_angle(value, "CB" , "CG" , "SD" , "C5'")
		beta  = dihedral_angle(value, "CG" , "SD" , "C5'", "C4'")
		gamma = dihedral_angle(value, "SD" , "C5'", "C4'", "O4'")
		delta = dihedral_angle(value, "C5'", "C4'", "C3'", "O3'")
		chi   = dihedral_angle(value, "O4'", "C1'", "N9" , "C2'")
		puts [pdb_id, key, alpha, beta, gamma, delta, chi].join(',')
	end
end

# main
puts ['pdb', 'chainID', 'alpha', 'beta', 'gamma', 'delta', 'chi'].join(',')

# filter file which ends with "pdb"
pdblist = Dir.entries("sam").delete_if { |element| element[-3,3] != 'pdb'}

pdblist.each do |f|
	calculate_sam("sam/"+f)
end