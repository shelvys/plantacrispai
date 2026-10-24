# Third-Party Licenses

PlantaCRISPAI is distributed under a proprietary license (see [`LICENSE`](LICENSE)).
It incorporates the third-party open source components listed below, each of
which remains governed by its own license. Nothing in PlantaCRISPAI's license
limits the rights granted to you under these licenses.

This file is required attribution. Do not remove it when redistributing any
build output.

## Summary

| Component | Version | License | Copyright |
|-----------|---------|---------|-----------|
| [react](https://github.com/facebook/react) | 19.x | MIT | Meta Platforms, Inc. and affiliates |
| [react-dom](https://github.com/facebook/react) | 19.x | MIT | Meta Platforms, Inc. and affiliates |
| [scheduler](https://github.com/facebook/react) | 0.28.x | MIT | Meta Platforms, Inc. and affiliates |
| [react-scripts](https://github.com/facebook/create-react-app) | 5.0.1 | MIT | Facebook, Inc. and its affiliates |
| [@babel/runtime](https://github.com/babel/babel) | 7.x | MIT | The Babel Team |

`react-scripts` is a build-time dependency. It is not shipped in the production
bundle, but it transitively installs a large dependency tree; run
`npm run licenses` (see below) for the complete, machine-generated inventory of
everything resolved in `node_modules`.

Only `react`, `react-dom`, `scheduler` and a small number of `@babel/runtime`
helpers are compiled into the deployed bundle. Their license notices are
preserved automatically by webpack in `build/static/js/*.license.txt`.

## Web fonts

The application loads **Inter** and **Fira Code** from Google Fonts at runtime.
Both are licensed under the [SIL Open Font License 1.1](https://openfontlicense.org/).
No font files are redistributed in this repository. The application degrades
gracefully to system fonts when the request is blocked or offline.

## Scientific methods

Published algorithms adapted by this software are credited in
[`docs/SCIENTIFIC_METHODS.md`](docs/SCIENTIFIC_METHODS.md). Those citations are
scientific attribution, not software licenses.

## Regenerating the full inventory

```bash
npx license-checker --production --summary
npx license-checker --json > third-party-inventory.json
```

---

## Full license texts

### MIT License — React, React DOM, Scheduler

MIT License

Copyright (c) Meta Platforms, Inc. and affiliates.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

### MIT License — react-scripts (Create React App)

MIT License

Copyright (c) 2013-present, Facebook, Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

### MIT License — @babel/runtime

MIT License

Copyright (c) 2014-present Sebastian McKenzie and other contributors

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
