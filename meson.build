project('cephe', 'c', license: 'GPL3+', version : '2.0')

#mdep = dependency('m', method: 'system')
mdep = meson.get_compiler('c').find_library('m', required : true)

executable('ephe',
  ['main.c', 'vect.c', 'body.c', 'orbit.c',
  'observ.c', 'instant.c', 'ephe.c'],
  dependencies: mdep,
  install: true)

