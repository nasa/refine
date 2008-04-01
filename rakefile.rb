
task :default => [:test]

task :test do
  system(' ( cd src && ./test_runner.rb ) ')
end
