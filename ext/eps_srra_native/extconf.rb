require 'mkmf'

RbConfig::MAKEFILE_CONFIG['CC'] = ENV['CC'] if ENV['CC'] # Note: this line is not strictly necessary, but gives the ability to easily use alternate compilers to build the extension. Since mkmf detects the compiler at Makefile creation time, this isn't very interesting until you consider static analysis tools, which tend to substitute the standard compiler with their own enhanced version. By having this line of code at the top, the Rakefile is prepared to let these static analysis tools do their thing (and help improve your code).

extension_name = 'eps_srra_native'

# - requiried development package -
# - here we do not have such pkg -
#unless pkg_config('library_to_link_to')
#  raise "library_to_link_to not found"
#end
#have_func('useful_function', 'library_to_link_to/lib.h')
#have_type('useful_type', 'library_to_link_to/lib.h')

create_header
create_makefile extension_name
