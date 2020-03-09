require 'simplecov'
require 'coveralls'

SimpleCov.formatter = Coveralls::SimpleCov::Formatter
SimpleCov.start do
   add_filter 'Examples'
   add_filter 'PubSrc'
   add_filter 'build'
   add_filter 'Doxumentation'
   add_filter 'Project/Parallel'
   add
end
