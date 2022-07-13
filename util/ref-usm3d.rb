#!/usr/bin/env ruby

class FUN3D

 @@dsl_methods = %w[ 
     fun3d_directory
     refine_directory
     usm3d_executable
     root_project
     number_of_processors
     mpirun_command
     iteration
     project_iteration_format
     first_iteration
     last_iteration
     window
     windows
     all_cl
     flo_cl
     adj_cl
     rad_cl
     ref_cl
     all_nl
     flo_nl
     adj_nl
     rad_nl
     pause
     breadcrumb_filename
     schedule_filename
     schedule_initial_complexity
     schedule_relative_tol
     schedule_complexity_growth
     schedule_subiteration_limit
     schedule_max_cores
     schedule_complexity_per_core
     schedule_max_complexity
     schedule_clean
 ]

 @@dsl_methods.each do |dsl_method|
  eval "def #{dsl_method}(value=@#{dsl_method});@#{dsl_method}=value;end"
 end

 DEFAULTS = <<-END_OF_DEFAULTS
  root_project ''
  project_iteration_format '%02d'
  first_iteration 1
  last_iteration 500
  pause 0
  schedule_initial_complexity 100000
  def schedule_force; cd; end
  schedule_relative_tol 0.001
  schedule_meshes_to_double 1
  schedule_complexity_growth 2.0
  schedule_subiteration_limit 20
  schedule_max_cores 400 # may change next line
  schedule_max_cores_inferred_from_batch_environment
  schedule_complexity_per_core 1000
  schedule_max_complexity 30_000_000
  schedule_clean true
  schedule_filename 'f3d_criteria.dat'
  breadcrumb_filename 'f3d.breadcrumbs'
  mpirun_command 'mpiexec'
 END_OF_DEFAULTS

 GROWABLE_DEFAULTS = <<-END_OF_GLOWABLE_DEFAULTS
  all_cl ""
  flo_cl ""
  adj_cl ""
  rad_cl ""
  ref_cl ""
  all_nl Hash.new
  flo_nl Hash.new
  adj_nl Hash.new
  rad_nl Hash.new
 END_OF_GLOWABLE_DEFAULTS

 def initialize()
  fun3d_directory File.dirname(__FILE__)
  refine_directory ''
  if (ENV::has_key?('USMEXE'))
    usm3d_executable ENV['USMEXE']
  else
    usm3d_executable 'usm3d.x'
  end
  eval DEFAULTS
  iteration 1
  eval GROWABLE_DEFAULTS
  load_case_specifics
 end

 def load_case_specifics(file_name = 'case_specifics')
  if File.exist?(file_name)
   eval GROWABLE_DEFAULTS
   eval IO.readlines(file_name).join("\n") 
  end
 end

 def project(iter=@iteration,win=@window)
  prj = @root_project
  prj += sprintf(project_iteration_format,iter) if iter
  prj += sprintf('_win'+project_iteration_format,win) if win
  prj
 end

 def input(hash = Hash.new, new_file = 'Flow/fun3d.nml')
  hash['project_rootname'] = "'#{project}'" unless hash['project_rootname']
  File.open(new_file,'w') do |f|
   IO.readlines('fun3d.nml').each do |line|
    hash.each_pair do |key,value|
     regex = /^\s*#{Regexp.escape(key)}\s*=.*/
     replacement = "#{key} = #{value}"
     line.gsub!(regex,replacement)
    end
    f.puts line
   end
  end
  `cp perturb.input Flow` if File.exist?( 'perturb.input' )
 end

 def have_project_extension?(extension)
   File.size?( "Flow/#{project}#{extension}" )
 end

 def have_file_extension?(suffix)
  if suffix =~ /adj/
    File.size?( "Adjoint/#{project}.adjoint" )
  else
    File.size?( "Flow/#{project}.#{suffix}" )
  end
 end

 def required_file_suffix(suffix)
   unless have_file_extension?(suffix)
     raise("#{project} #{suffix} file zero length or missing")
   end
 end

 def expect_file(filename)
   unless File.size?(filename)
     raise("expected file #{filename} zero length or missing")
   end
 end

 def exec_prefix()
   cores = number_of_processors ||
           schedule_cores(schedule_max_cores,schedule_complexity_per_core)
   if (cores > 1) then
     "#{mpirun_command} -np #{cores}"
   else
     " "
   end
 end

 def schedule_current()
   complexity = schedule_initial_complexity
   if File.file?(schedule_filename)
     lines = IO.readlines(schedule_filename)
     complexity = lines.last.split[1].to_i # last complexity
   end
   return complexity
 end
 
 def schedule_parse()
   current_complexity = schedule_initial_complexity
   next_complexity = schedule_initial_complexity
   if File.file?(schedule_filename)
     lines = IO.readlines(schedule_filename)
     current_complexity = lines.last.split[1].to_i # last line complexity
     next_complexity = lines.last.split[2].to_i # last line complexity target
   end
   return current_complexity, next_complexity
 end

 def first_iteration_inferred_from_newest_meshb
   if (Dir['Flow/*.meshb'].empty?) then
     first_iteration 1
   else
     newest_meshb = Dir['Flow/*.meshb'].sort_by{|meshb| File.mtime(meshb)}.last
     iter = newest_meshb.sub(/Flow\/#{root_project}/,'').sub(/\.meshb/,'').to_i
     first_iteration iter
   end
 end

 def schedule_max_cores_inferred_from_batch_environment
   if (ENV::has_key?('PBS_NODEFILE'))
     schedule_max_cores IO.readlines(ENV['PBS_NODEFILE']).length
   end
   if (ENV::has_key?('SLURM_JOB_NUM_NODES'))
     schedule_max_cores ENV['SLURM_JOB_NUM_NODES'].to_i
   end
 end

 def schedule_cores(max_cores = schedule_max_cores,
                    complexity_per_core = schedule_complexity_per_core)
   cores_desired =
     (schedule_current().to_f / schedule_complexity_per_core.to_f).round
   return [2,[cores_desired,schedule_max_cores].min].max
 end

 def schedule_append(current_complexity, next_complexity)
   open(schedule_filename, 'a') do |f|
     f.puts "#{iteration} #{current_complexity} #{next_complexity} #{schedule_force}"
   end
 end

 def schedule_parse_history(current)
   iterations_at_complexity = 0
   forces = []
   if File.file?(schedule_filename)
     lines = IO.readlines(schedule_filename)
     lines.reverse.each do |line|
       complexity = line.split[1].to_i
       if (complexity != current) then
         break
       end
       iterations_at_complexity += 1
       forces << line.split[3].to_f if (iterations_at_complexity <= 2)
     end
   end
   iterations_at_complexity += 1
   forces << schedule_force()
   return iterations_at_complexity, forces
 end

 def schedule_next()
  current_complexity, next_complexity = schedule_parse()
  iterations_at_complexity, forces = schedule_parse_history(current_complexity)
  puts " #{iterations_at_complexity} meshes at #{schedule_current()}"
  puts " force history #{forces.join(' ')}"
  if (iterations_at_complexity >= schedule_subiteration_limit) then
    next_complexity *= 2
  else
    if (iterations_at_complexity > 3) then
      force_range = forces.max - forces.min
      force_mid = 0.5*(forces.max + forces.min)
      relative_tol = schedule_relative_tol
      printf("mid %f range %e tol %e\n",
             force_mid, force_range, relative_tol*force_mid)
      next_complexity *= 2 if ( relative_tol*force_mid.abs > force_range)
    end
  end
  current_complexity=[next_complexity,
    (current_complexity.to_f*schedule_complexity_growth.to_f).ceil].min
  schedule_append(current_complexity, next_complexity)
  return current_complexity
 end

 def schedule_meshes_to_double( number_of_meshes )
   schedule_complexity_growth 2.0**(1.0/number_of_meshes.to_f)
 end
 
 def nap( duration = @pause )
   if duration > 0 
     puts "sleeping #{duration} seconds"
     sleep(duration)
   end
 end

 def sh(comm)
   puts comm
   nap
   start_time = Time.now
   raise("nonzero exit code") unless system(comm)
   printf("%d %f # sh\n",@iteration,Time.now-start_time)
 end

 def viscous_tags( mapbc = "Flow/#{project}.mapbc", separator = ',' )
   wall_distance_varietals = [4000,-4000,4075,4100,4110,-4110,-4100,6200,6210]
   tags = Array.new
   lines = IO.readlines(mapbc)
   n = lines[0].to_i
   n.times.each do |i|
     cols = lines[1+i].split()
     face_tag = cols[0].to_i
     bc_type = cols[1].to_i
     tags << face_tag if ( wall_distance_varietals.include?(bc_type) )
   end
   tags.join(separator)
 end

 def dist(command = nil, stdoutfile = 'dist_out', expected_file_without_path = nil)
  `cp #{root_project}.mapbc Flow/#{project}.mapbc` if File.exist?("#{root_project}.mapbc")
  `cp #{project(1)}.mapbc Flow/#{project}.mapbc` if File.exist?("#{project(1)}.mapbc")
   input_file = "#{project}.meshb"
   input_file = "#{project}.b8.ugrid" if (File.exist?("Flow/#{project}.b8.ugrid"))
   input_file = "#{project}.lb8.ugrid" if (File.exist?("Flow/#{project}.lb8.ugrid"))
   capture = " < /dev/null | tee #{stdoutfile} > #{project}_#{stdoutfile}"
   comm = command || "( cd Flow && " +
             exec_prefix +
             " ParallelDistanceCalculator " +
             " #{input_file} " +
             " --commas #{viscous_tags()} " +
             capture +
             " )"
   puts comm
   nap
   raise("ParallelDistanceCalculator nonzero exit code") unless system(comm)
   nap
   expected_file_without_path =
     expected_file_without_path || #{project}-distance.solb"
   expect_file("Flow/#{expected_file_without_path}")
 end

 def clean_step()
   clean_files = []
   %w[ .lb8.ugrid -distance.solb -distance.solb.names .flow .grid_info .mapbc .meshb -restart.solb _volume.solb].each do |suffix|
     clean_files << "Flow/#{project}#{suffix}"
   end
   sh("rm -rf #{clean_files.join(' ')}")
 end

 def ref_loop_or_exit(opts = '')
   current_complexity = schedule_current()
   complexity = schedule_next()
   complexity_limit = schedule_max_complexity
   if (complexity > complexity_limit) then
     puts "next complexity #{complexity} larger than complexity_limit #{complexity_limit}"
     puts "done."
     exit
   end
   ref("refmpi loop #{project} #{project(@iteration+1)} #{complexity} "+
       ref_cl+' '+opts,
       "refine_out",
       "#{project(@iteration+1)}.meshb")
   if (schedule_clean && current_complexity == complexity) then
     clean_step()
   end
 end

 def refdist(opts = '')
   `cp #{root_project}.mapbc Flow/#{project}.mapbc` if File.exist?("#{root_project}.mapbc")
   `cp #{project(1)}.mapbc Flow/#{project}.mapbc` if File.exist?("#{project(1)}.mapbc")
   input_file = "#{project}.meshb"
   input_file = "#{project}.b8.ugrid" if (File.exist?("Flow/#{project}.b8.ugrid"))
   input_file = "#{project}.lb8.ugrid" if (File.exist?("Flow/#{project}.lb8.ugrid"))
   ref("refmpi distance " +
       " #{input_file} " +
       " #{project}-distance.solb " +
       " --fun3d-mapbc #{project}.mapbc ",
       "refdist_out",
       "#{project}-distance.solb")
 end

 def ref( comm, stdoutfile = '', expected_file_without_path = nil)
   capture = (stdoutfile.empty? ? '':" < /dev/null | tee #{stdoutfile} > #{project}_#{stdoutfile}")
   command = "( cd Flow && " +
             exec_prefix +
             " " +
             @refine_directory +
             comm +
             capture +
             " )"
   puts command
   nap
   start_time = Time.now
   raise("refine nonzero exit code") unless system(command)
   printf("%d %f # ref\n",@iteration,Time.now-start_time)
   nap
   if ( expected_file_without_path ) then
     expect_file('Flow/'+expected_file_without_path)
   else
     last_arg_without_leading_dash =
       comm.split(' ').delete_if{ |f| f.match(/^-/) }.last
     expect_file('Flow/'+last_arg_without_leading_dash)
   end
 end

 def flo( extra_nl = Hash.new, extra_cl = '' )
  `cp #{project(1)}.knife Flow/#{project}.knife` if File.exist?("#{project(1)}.knife")
  `cp #{project(1)}.cutbc Flow/#{project}.cutbc` if File.exist?("#{project(1)}.cutbc")
  `cp sfe.cfg Flow` if File.exist?("sfe.cfg")
  name_list = Hash.new
  name_list['restart_read'] = "'#{have_file_extension?('flow') ? 'on':'off'}'"
  name_list['distance_from_file'] = "'#{have_project_extension?('-distance.solb') ? "#{project}-distance.solb":''}'"
  name_list['import_from'] = "'#{have_project_extension?('-restart.solb') ? "#{project}-restart.solb":''}'"
  name_list = name_list.merge(all_nl).merge(flo_nl).merge(extra_nl)
  input(name_list)
  system("( cd Flow && cp fun3d.nml #{project}_flow_fun3d.nml)")
  system("( cd Flow && cp sfe.cfg #{project}_flow_sfe.cfg)") if File.exist?("Flow/sfe.cfg")
  `cp presb.input Flow` if File.exist?("presb.input")
  `cp SBtarget.sig Flow` if File.exist?("SBtarget.sig")
  `cp #{root_project}.mapbc Flow/#{project}.mapbc` if File.exist?("#{root_project}.mapbc")
  `cp #{project(1)}.mapbc Flow/#{project}.mapbc` if File.exist?("#{project(1)}.mapbc")
  command = "( cd Flow && " +
            exec_prefix +
            " nodet_mpi " +
             [all_cl,flo_cl,extra_cl].join(' ') +
            "< /dev/null | tee flow_out > #{project}_flow_out )"
  puts command
  nap
  start_time = Time.now
  raise("flow solver nonzero exit code") unless system(command)
  printf("%d %f # flo\n",@iteration,Time.now-start_time)
  nap
  read_forces( )
 end

 def complex( )
  `mkdir Complex`
  %w[ fgrid fastbc cogsg bc mapbc knife cutbc ugrid msh amdba ].each do |ext| 
   file = "#{project}.#{ext}"
   if File.exist?(file)
    command = "cp #{file} Complex"
    puts command
    system command
   end
  end
  Dir["Flow/fun3d.nml"].each do |target|
   `( cd Complex && ln -s ../#{target} .)`
  end
  `cp perturb.input Complex`
  command = "( cd Complex && " +
            exec_prefix +
            " complex_nodet_mpi " +
             [all_cl,flo_cl].join(' ') +
            "< /dev/null | tee flow_out > #{project}_flow_out )"
  puts command
  nap
  start_time = Time.now
  raise("complex flow solver nonzero exit code") unless system(command)
  printf("%d %f # complex\n",@iteration,Time.now-start_time)
  nap
 end

 def read_forces( )
  IO.readlines("Flow/#{project}.forces").each do |line|
   line.scan(/C\S+\s*=\s*\S+/) do |force|
    eval("def #{force.gsub(/=/,';')};end".downcase)
   end
  end
 end

 def adj( extra_nl = Hash.new, extra_cl = '' )
  required_file_suffix('flow')
  name_list = Hash.new.merge(all_nl).merge(adj_nl).merge(extra_nl)
  name_list['restart_read'] = "'#{ have_file_extension?('adj') ? 'on':'off'}'"
  input(name_list)
  system("( cd Flow && cp fun3d.nml #{project}_dual_fun3d.nml)")
  `cp presb.input Adjoint` if File.exist?("presb.input")
  `cp SBtarget.sig Adjoint` if File.exist?("SBtarget.sig")
  command = "( cd Adjoint && " +
            exec_prefix +
            " dual_mpi " +
            [all_cl,adj_cl,extra_cl].join(' ') +
            "< /dev/null | tee dual_out > #{project}_dual_out " +
            ")"
  puts command
  nap
  start_time = Time.now
  raise("adjoint solver nonzero exit code") unless system(command)
  printf("%d %f # adj\n",@iteration,Time.now-start_time)
  nap
  %w[ pressure_signatures.dat sboom.data SBground.sig Gradient.plt ].each do |file|
   `cp Adjoint/#{file} Adjoint/#{project}_#{file}` if File.exist?( "Adjoint/#{file}" )
  end
 end

 def rad( extra_nl = Hash.new, extra_cl = '' )
  required_file_suffix('flow')
  required_file_suffix('adj')
  name_list = Hash.new.merge(all_nl).merge(rad_nl).merge(extra_nl)
  name_list['adapt_project'] = "'#{project(@iteration+1)}'"
  input(name_list)
  system("( cd Flow && cp fun3d.nml #{project}_rad_fun3d.nml)")
  `cp faux_input Adjoint` if File.exist?("faux_input")
  `cp #{project(1)}.freeze Adjoint/#{project}.freeze` if File.exist?("#{project(1)}.freeze")
  `cp #{project(1)}.knife Adjoint/refine.knife` if File.exist?("#{project(1)}.knife")
  `cp presb.input Adjoint` if File.exist?("presb.input")
  `cp SBtarget.sig Adjoint` if File.exist?("SBtarget.sig")
  command = "( cd Adjoint && " +
            exec_prefix +
            " dual_mpi " +
            " --rad " +
            " --adapt " +
            [all_cl,rad_cl,extra_cl].join(' ')  +
            "< /dev/null | tee rad_out > #{project}_rad_out " +
            ")"
  puts command
  nap
  start_time = Time.now
  raise("rad adapt nonzero exit code") unless system(command)
  printf("%d %f # rad\n",@iteration,Time.now-start_time)
  nap
 end

 def adapt( extra_nl = Hash.new, extra_cl = '' )
  required_file_suffix('flow')
  name_list = Hash.new
  name_list['adapt_project'] = "'#{project(@iteration+1)}'"
  name_list = name_list.merge(all_nl).merge(rad_nl).merge(extra_nl)
  name_list['restart_read'] = "'on'"
  input(name_list)
  system("( cd Flow && cp fun3d.nml #{project}_adapt_fun3d.nml)")
  `cp faux_input Flow` if File.exist?("faux_input")
  `cp #{project(1)}.freeze Flow/#{project}.freeze` if File.exist?("#{project(1)}.freeze")
  command = "( cd Flow && " +
            exec_prefix +
            " nodet_mpi " +
            " --adapt " +
            [all_cl,rad_cl,extra_cl].join(' ')  +
            "< /dev/null | tee adapt_out > #{project}_adapt_out " +
            ")"
  puts command
  nap
  start_time = Time.now
  raise("adapt nonzero exit code") unless system(command)
  printf("%d %f # adapt\n",@iteration,Time.now-start_time)
  nap
 end

 def usm( )
  sh("cp #{root_project}.inpt Flow/#{project}.inpt") if File.exist?("#{root_project}.inpt")
  sh("cp #{project(1)}.inpt Flow/#{project}.inpt") if File.exist?("#{project(1)}.inpt")
  sh("cp #{root_project}.mapbc Flow/#{project}.mapbc") if File.exist?("#{root_project}.mapbc")
  sh("cp #{project(1)}.mapbc Flow/#{project}.mapbc") if File.exist?("#{project(1)}.mapbc")

  sh("rm -rf Flow/tet.out Flow/timer.out Flow/fort.* Flow/*.finaliface.* " +
     " Flow/*.frincells.* Flow/*.fringebdry.*  Flow/*.parpart.* " +
     " Flow/*.walldist.* Flow/CELLS*")
  sh("rm -rf Flow/volume.plt Flow/surface.plt Flow/hist.plt")

  command = "( cd Flow && " +
            exec_prefix +
            " #{usm3d_executable} #{project} " +
            "< /dev/null | tee flow_out > #{project}_flow_out )"
  puts command
  nap
  start_time = Time.now
  raise("flow solver nonzero exit code") unless system(command)
  printf("%d %f # flo\n",@iteration,Time.now-start_time)
  nap

  expect_file("Flow/volume.plt")
  sh("mv Flow/volume.plt Flow/#{project}_volume.plt")
  sh("mv Flow/surface.plt Flow/#{project}_surface.plt")
  sh("mv Flow/hist.plt Flow/#{project}_hist.plt")

  # read forces for scheduling complexity
  line = IO.readlines("Flow/#{project}_hist.plt").last
  col = line.split(' ')
  eval("def res;    #{col[2]}; end")
  eval("def logres; #{col[3]}; end")
  eval("def cl;     #{col[4]}; end")
  eval("def cd;     #{col[5]}; end")
  eval("def cdv;    #{col[6]}; end")
  eval("def cm;     #{col[7]}; end")
 end

 def global_refine
   flo(extra_nl={}, extra_cl="--embed_grid")
   %w[ fgrid fastbc cogsg bc mapbc knife cutbc ugrid lb8.ugrid b8.ugrid msh amdba ].each do |ext| 
     file = "#{project}_embed.#{ext}"
     if File.exist?("Flow/#{file}")
       command = "( cd Flow && mv #{file} #{project(@iteration+1)}.#{ext} )"
       puts command
       system command
     end
   end
 end

 def iteration_steps
   refdist
   flo
   ref_loop_or_exit
 end

 def iterate
   @iteration = first_iteration
   setup if ( 1 == @iteration )
   while ( @iteration <= last_iteration )
     load_case_specifics
     iteration_steps
     break if ( @iteration >= last_iteration )
     @iteration += 1
   end
 end

 def full_path_for(exec)
   ['','../FUN3D_90','../Adjoint','../Complex/FUN3D_90','../utils'].each do |dir|
     path = File.join(File.join(fun3d_directory,dir),exec)
     return path if ( File.executable? path )
   end
  puts " Unable to find #{exec} in path."
  puts " Make sure #{$0} is in the FUN3D installed bin directory,"
  puts "   not a copy or link."
  puts " Alternatively, set fun3d_directory in case_specifics."
  raise("can not find #{exec} in path")
 end

 def setup
  `mkdir -p Flow Adjoint`
  %w[ fgrid fastbc cogsg bc mapbc knife cutbc meshb ugrid lb8.ugrid b8.ugrid msh amdba ].each do |ext|
   file = "#{project}.#{ext}"
   if File.exist?(file)
    command = "cp #{file} Flow"
    puts command
    system command
   end
  end

  gengas_files = ["tdata", "species_thermo_data", "species_transp_data",
                  "species_transp_data_0", "kinetic_data"]
  optional_files = ["#{root_project}_g.msh", "moving_body.input",
                    "rotor.input"] + gengas_files
  optional_files.each do |file|
   if File.exist?(file)
    command = "cp #{file} Flow"
    puts command
    system command
   end
  end

  100.times do |prop|
   file = sprintf("propeller_properties%d.dat",prop)
   if File.exist?(file)
    command = "cp #{file} Flow"
    puts command
    system command
   end
  end

 end

end

def conditional_tail(file, lines=5)
  if File.exists?(file)
    puts file+' -------------'
    system("tail -#{lines} #{file}")
  end
end

def rubber_data_writer(function)
  File.open('rubber.data','w') do |f|
    f.puts <<END_OF_RUBBER
################################################################################
########################### Design Variable Information ########################
################################################################################
Global design variables (Mach number, AOA, Yaw, Noninertial rates)
  Var Active         Value               Lower Bound            Upper Bound
 Mach    0   0.100000000000000E+01  0.000000000000000E+00  0.500000000000000E+01
  AOA    0   0.100000000000000E+01  0.000000000000000E+00  0.500000000000000E+01
  Yaw    0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
xrate    0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
yrate    0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
zrate    0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
Number of bodies
    0
#########################################################################
############################ Function Information #######################
#########################################################################
Number of composite functions for design problem statement
    1
#########################################################################
Cost function (1) or constraint (2)
    1
If constraint, lower and upper bounds
              0.100000000000000              0.500000000000000
Number of components for function   1
    1
Physical timestep interval where function is defined
    1     1
Composite function weight, target, and power
1.0 0.0 1.0
Components of function   1: boundary id (0=all)/name/value/weight/target/power
    0 #{function}                     -0.0           1.000    0.00000  1.00
Current value of function   1
            -0.000000000000000
Current derivatives of function wrt global design variables
             0.000000000000000
             0.000000000000000
             0.000000000000000
             0.000000000000000
             0.000000000000000
             0.000000000000000
END_OF_RUBBER
  end
end

if ( __FILE__ == $0 )

  STDOUT.sync = true # to encourage real-time status messages

  case ARGV[0]
  when 'start'
    system("nohup #{$0} iterate < /dev/null >> output &")
    puts 'started, see `output\' file or f3d watch'
  when 'view'
    system('pwd')
    conditional_tail 'output', 3
    conditional_tail 'Flow/flow_out'
    conditional_tail 'Flow/refine_out'
    conditional_tail 'Flow/adapt_out'
    conditional_tail 'Adjoint/dual_out'
    conditional_tail 'Adjoint/rad_out'
  when 'watch'
    exec("watch #{$0} view")
  when 'hist'
    exec(File.join(File.dirname(__FILE__),'hist2gnuplot'))
  when 'clean'
    fun3d = FUN3D.new
    fun3d.sh("rm -rf Flow Adjoint Complex output "+
             "#{fun3d.breadcrumb_filename} #{fun3d.schedule_filename}")
  when 'shutdown'
    comm = 'killall -9 nodet_mpi dual_mpi complex_nodet_mpi ruby' ; puts comm ; system comm
  when 'check'
    fun3d = FUN3D.new
    fun3d.setup
    fun3d.flo
    fun3d.adj
    fun3d.complex
    comm = 'grep dfunc Complex/flow_out Adjoint/dual_out'
    puts comm ; system comm
  when 'iterate'
    fun3d = FUN3D.new
    File.open(fun3d.breadcrumb_filename,'w') do |f|
      f.puts "#{`uname -n`.chomp}:#{Process.pid}"
    end
    fun3d.iterate
  when 'function'
    rubber_data_writer(ARGV[1]||'cd')
  else
    puts "usage: #{File.basename($0)} <command>"
    puts
    puts " <command>       description"
    puts " ---------       -----------"
    puts " start           Start adaptation (in background)"
    puts " iterate         Start adaptation (blocking)"
    puts " view            Echo a single snapshot of stdout"
    puts " watch           Watch the result of view"
    puts " shutdown        Kill all running fun3d and ruby processes"
    puts " clean           Remove output and sub directories"
    puts " function [name] write rubber.data with cost function [name]"
    puts
    puts " Defaults in the case_specific input file:"
    puts FUN3D::DEFAULTS
    puts FUN3D::GROWABLE_DEFAULTS
  end
end

