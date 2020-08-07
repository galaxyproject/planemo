<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>{{ title }}</title>

    <!-- Bootstrap -->
    {{ bootstrap_style }}
    {{ custom_style }}

  </head>
  <body>

    <nav class="navbar navbar-inverse navbar-fixed-top" role="navigation">
      <div class="container-fluid">
        <div class="navbar-header">

          <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false" aria-controls="navbar">
            <span class="sr-only">Toggle navigation</span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
          </button>
          <a class="navbar-brand" href="#">{{ title }}</a>
        </div>
        <div id="navbar" class="navbar-collapse collapse">
          <ul class="nav navbar-nav navbar-right">
            <li><a href="https://galaxyproject.org">Galaxy</a></li>
            <li><a href="https://planemo.readthedocs.org">Planemo</a></li>
          </ul>
          <div class="navbar-form navbar-right">
          </div>
        </div>
      </div>
    </nav>

    <div class="container-fluid">
      <div id="overview-content" class="row col-md-offset-1">
      </div>
    </div>

    {{ jquery_script }}
    {{ bootstrap_script }}
    <style>
details > summary {
  vertical-align: text-top;
}
     </style>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/showdown/1.9.1/showdown.min.js" integrity="sha512-L03kznCrNOfVxOUovR6ESfCz9Gfny7gihUX/huVbQB9zjODtYpxaVtIaAkpetoiyV2eqWbvxMH9fiSv5enX7bw==" crossorigin="anonymous"></script>
    <script>
        showdown.setFlavor('github');
        showdown.setOption('emoji', true);
        var converter = new showdown.Converter();
        var test_data = '{{ raw_data }}';
        var html = converter.makeHtml(atob(test_data));
        var target = document.getElementById('overview-content');
        target.innerHTML = html;
    </script>
  </body>
</html>
