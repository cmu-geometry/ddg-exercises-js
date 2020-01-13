## Building the Documentation

This documentation was generated using [jsdoc](http://usejsdoc.org/) and [ink-docstrap](https://www.npmjs.com/package/ink-docstrap). We used a modified version of the [cosmo](http://docstrap.github.io/docstrap/themes/cosmo/) theme. To build the documentation, first instal jsdoc and ink-docstrap with npm in the parent directory of `ddg-exercises-js`. Then copy `doc-config/site.cosmo-rohan.css` into `node_modules/ink-docstrap/template/static/styles/`. Now you can build the documentation by running
```
jsdoc -c ddg-exercises-js/doc-config/jsdoc.conf.json -t node_modules/ink-docstrap/template/ -R ddg-exercises-js/README.md -r -d ddg-exercises-js/docs ddg-exercises-js/
```
