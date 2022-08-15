project('cephe', 'c', license: 'GPL3+', version : '2.0')

#mdep = dependency('m', method: 'system')
mdep = meson.get_compiler('c').find_library('m', required : true)

intldep = dependency('intl', required: false)

conf_data = configuration_data()
conf_data.set('HAVE_LIBINTL', intldep.found())
configure_file(output : 'config.h', configuration: conf_data)
confinc = include_directories('.')

executable('ephe',
  ['main.c', 'vect.c', 'body.c', 'orbit.c',
  'observ.c', 'instant.c', 'ephe.c'],
  dependencies: [mdep, intldep],
  include_directories: confinc,
  install: true)

